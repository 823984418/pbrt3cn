# 层次包围盒结构

层次包围盒结构（BVH）是一种基于图元划分的光线求交加速的方法，其中图元被划分为不相交的层次结构。（相比之下，空间划分通常将空间划分为不相交的层次结构。图4.2展示了简单场景的层次包围盒结构。图元存储在树叶中，每个节点储存了其下方的节点中图元的包围盒。因此，当光线穿过树时，只要它不与节点的包围盒相交，都可以跳过该节点下方的子树。

<img src="https://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Primitives and hierarchy.svg">图 4.2：简单场景的层次包围盒结构。（a） 图元的少量集合，用虚线显示边界框。图元根据邻近性进行聚合;在这里，球体和等边三角形由另一个边界框边界，然后由包含整个场景的（均以实线显示）包围盒。（b） 相应的层次包围盒结构。根节点保留整个场景的包围盒。在这里，它有两个子项，一个存储一个包含球体和等边三角形的包围盒（这反过来又有这些图元作为其子对象），另一个存储包含瘦三角形的包围盒。</img>

图元划分的一个特性是每个图元仅在层次结构中出现一次。相反，空间划分中图元可能与空间重叠于多个区域，因此当光线穿过它们时，需要多次测试交点。此属性的另一个含义是表示图元划分层次结构所需的内存量是有边界的。对于在每个叶中存储单个基元的二进制层次包围盒结构，节点总数为$2n - 1$其中$n$是图元数。将有$n$个叶子节点和$n - 1$个枝节点。如果叶子节点存储多个图元，则需要更少的节点。

层次包围盒的构建效率高于K维空间划分树，K维空间划分树通常提供比层次包围盒稍快的射线求交测试，但构建时间要长得多。另一方面，层次包围盒在数值上通常更健壮，而且比K维空间划分树更不容易因舍入错误而错过交点。

BVH加速器，`BVHAccel`被定义于`accelerators/bvh.h`和`accelerators/bvh.cpp`。除了要存储的基元和任何叶节点中可以存储的最大基元数之外，其构造函数还采用一个枚举值，该值描述了在划分图元以生成树时要使用的四种算法中的哪一种。默认值`SAH`表示应使用第4.3.2节中讨论的基于"表面积启发式"的算法。另一种选择，`HLBVH`，这是在第4.3.3节中讨论的，可以更有效地构造（更容易并行化），但它不能构建像SAH那样有效的树。其余两种方法使用更少的计算来生成树，但创建质量相当低的树。

```
<<BVHAccel公共类型>>=
enum class SplitMethod {
    SAH, HLBVH, Middle, EqualCounts
};

<<BVHAccel方法定义>>=
BVHAccel::BVHAccel(const std::vector<std::shared_ptr<Primitive>> &p,
         int maxPrimsInNode, SplitMethod splitMethod)
     : maxPrimsInNode(std::min(255, maxPrimsInNode)), primitives(p),
       splitMethod(splitMethod) {
    if (primitives.size() == 0)
        return;
    <<从图元构建BVH>>=
       <<从所有图元中构建图元信息数组>>=
          std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
          for (size_t i = 0; i < primitives.size(); ++i)
              primitiveInfo[i] = { i, primitives[i]->WorldBound() };

       <<从使用图元信息数组的所有图元构建BVH树>>=
          MemoryArena arena(1024 * 1024);
          int totalNodes = 0;
          std::vector<std::shared_ptr<Primitive>> orderedPrims;
          BVHBuildNode *root;
          if (splitMethod == SplitMethod::HLBVH)
              root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
          else
              root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
                                    &totalNodes, orderedPrims);
          primitives.swap(orderedPrims);
          

       <<计算BVH树的深度优先遍历的表示>>=
          nodes = AllocAligned<LinearBVHNode>(totalNodes);
          int offset = 0;
          flattenBVHTree(root, &offset);
          
}

<<BVHAccel私有数据>>=
const int maxPrimsInNode;
const SplitMethod splitMethod;
std::vector<std::shared_ptr<Primitive>> primitives;

```

## BVH构造
在此实施中，BVH 构建分为三个阶段。首先，计算每个基元的边界信息并存储在将在树构造期间使用的数组中。接下来，使用在拆分方法中编码的算法选择构建树。结果是一个二叉树，其中每个内部节点都保存指向其子节点的指针，每个叶节点保存对一个或多个基元的引用。最后，此树将转换为更紧凑（从而更高效）的无指针表示，用于渲染期间。（此方法的实现更为简单，而不是直接在树构造期间计算无指针表示形式，这也是可能的。

```
<<从图元构建BVH>>= 
    <<从所有图元中构建图元信息数组>> 
    <<从使用图元信息数组的所有图元构建BVH树>> 
    <<计算BVH树的深度优先遍历的表示>> 
```

对于要存储在BVH中的每个图元，我们将其包围盒的中心、其完整包围盒及其索引存储在`BVHPrimitiveInfo`结构实例中的图元数组中。

```
<从所有图元中构建图元信息数组>>=
std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
for (size_t i = 0; i < primitives.size(); ++i)
    primitiveInfo[i] = { i, primitives[i]->WorldBound() };

<<BVHAccel本地声明>>= 
struct BVHPrimitiveInfo {
    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
        : primitiveNumber(primitiveNumber), bounds(bounds),
          centroid(.5f * bounds.pMin + .5f * bounds.pMax) { }
    size_t primitiveNumber;
    Bounds3f bounds;
    Point3f centroid;
};
```

层次结构构造现在可以开始。如果选择了`HLBVH`构造算法，则调用`HLBVHBuild()`来生成树。其他三种构造算法均由`recursiveBuild()`处理。这些函数的初始调用将传递要存储在树中的所有图元。它们返回指向树根的指针，该根目录用`BVHBuildNode`结构表示。树节点应与提供的`MemoryArena`一起分配，并且创建的总数应存储在`*totalNodes`中。

树构造过程的一个重要副作用是，通过`orderedPrims`参数返回指向基元的新的指针数组;此数组存储排序的基元，以便叶节点中的基元在数组中占据连续范围。在树构造后，它与原始基元数组交换。

```
<<从使用图元信息数组的所有图元构建BVH树>>= 
MemoryArena arena(1024 * 1024);
int totalNodes = 0;
std::vector<std::shared_ptr<Primitive>> orderedPrims;
BVHBuildNode *root;
if (splitMethod == SplitMethod::HLBVH)
    root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
else
    root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
                          &totalNodes, orderedPrims);
primitives.swap(orderedPrims);
```

