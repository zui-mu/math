# 向量空间、列空间与零空间

## 课程主线与当前范围

这部分课程始终围绕同一个线性方程展开：

```math
A\mathbf{x}=\mathbf{b},
\qquad
A:\mathbb{R}^n\longrightarrow\mathbb{R}^m.
```

其中 $`\mathbf{x}`$ 是输入，$`A\mathbf{x}`$ 是输出，$`\mathbf{b}`$ 是希望得到的目标输出。当前课程主要从两个方向理解矩阵 $`A`$：

| 观察方向 | 已学空间 | 描述的内容 |
| --- | --- | --- |
| 输出端 $`\mathbb{R}^m`$ | 列空间 $`\mathrm{Col}(A)`$ | $`A`$ 能产生哪些输出 |
| 输入端 $`\mathbb{R}^n`$ | 零空间 $`N(A)`$ | 哪些输入会被 $`A`$ 变成零向量 |

后文先建立子空间与张成，再从矩阵的列进入列空间；接着用消元识别主元、自由变量和秩，由此得到零空间，最后把列空间与零空间合起来描述 $`A\mathbf{x}=\mathbf{b}`$ 的全部解。

课程中已经计算并比较过 $`N(A^T)`$，但行空间、左零空间以及四个基本子空间之间的正交关系还没有正式展开。

## 一、向量空间与子空间

### 1. $`\mathbb{R}^n`$ 的含义

```math
\mathbb{R}^n
= \{(x_1,x_2,\ldots,x_n)\mid x_1,x_2,\ldots,x_n\in\mathbb{R}\}
```

$`\mathbb{R}^n`$ 表示所有由 $`n`$ 个实数组成的向量。比如，$`\mathbb{R}^2`$ 是所有二维实向量的集合，$`\mathbb{R}^3`$ 是所有三维实向量的集合。

### 2. 子空间与封闭性

子空间是原向量空间中的一个子集，并且它自己也构成向量空间。判断非空集合 $`S`$ 是否为子空间，可以检查：对任意 $`\mathbf{u},\mathbf{v}\in S`$ 和任意实数 $`a,b`$，是否都有

```math
a\mathbf{u}+b\mathbf{v}\in S.
```

这叫作对线性组合封闭。它同时包含两个基本要求：

```math
\mathbf{u}+\mathbf{v}\in S,
\qquad
c\mathbf{u}\in S.
```

令 $`a=b=0`$，还可以得到零向量，因此**子空间一定包含零向量**。

例如，$`x`$ 轴

```math
S=\{(x,0)\mid x\in\mathbb{R}\}
```

满足

```math
a(x_1,0)+b(x_2,0)=(ax_1+bx_2,0)\in S,
```

所以它是 $`\mathbb{R}^2`$ 的子空间。

水平直线

```math
T=\{(x,1)\mid x\in\mathbb{R}\}
```

不是子空间，因为它不包含 $`(0,0)`$，并且

```math
(2,1)+(5,1)=(7,2)\notin T.
```

封闭性始终针对当前正在判断的集合。把集合扩大后，新集合可能封闭，但这不能改变原集合不是子空间的事实。

### 3. 张成

一组向量所有线性组合构成的集合叫作它们的**张成**（span）：

```math
\mathrm{span}\{\mathbf{v}_1,\ldots,\mathbf{v}_k\}
= \{c_1\mathbf{v}_1+\cdots+c_k\mathbf{v}_k
\mid c_1,\ldots,c_k\in\mathbb{R}\}.
```

张成集合包含零向量，并且对加法和数乘封闭，所以它一定是子空间。几何上，它可能是原点、过原点的直线、过原点的平面，或者更高维的子空间。

子空间与张成为后面两个核心对象提供了共同语言：

**列空间是矩阵各列的张成，零空间是齐次方程全部解组成的子空间。**

## 二、从矩阵方程到列空间

### 1. $`A\mathbf{x}=\mathbf{b}`$ 的尺寸

若

```math
A\in\mathbb{R}^{m\times n},
```

那么

```math
\underbrace{A}_{m\times n}
\underbrace{\mathbf{x}}_{n\times1}
= \underbrace{\mathbf{b}}_{m\times1}.
```

矩阵结构与方程组的对应关系是：

| 矩阵结构 | 方程组含义 |
| --- | --- |
| $`m`$ 行 | $`m`$ 个方程 |
| $`n`$ 列 | $`n`$ 个未知数 |
| 第 $`i`$ 行 | 第 $`i`$ 个方程的系数 |
| 第 $`j`$ 列 | 未知数 $`x_j`$ 在全部方程中的系数 |

同一个矩阵可以从三个角度理解：

1. 按行看，它是一组方程。
2. 按列看，$`A\mathbf{x}`$ 是矩阵各列的线性组合。
3. 看成线性变换，它把 $`\mathbb{R}^n`$ 中的输入送到 $`\mathbb{R}^m`$ 中：

```math
A:\mathbb{R}^n\longrightarrow\mathbb{R}^m.
```

### 2. 矩阵乘向量的列视角

把矩阵写成列向量形式：

```math
A=
\begin{bmatrix}
\mathbf{a}_1&\mathbf{a}_2&\cdots&\mathbf{a}_n
\end{bmatrix}.
```

那么

```math
A\mathbf{x}
= x_1\mathbf{a}_1+x_2\mathbf{a}_2+\cdots+x_n\mathbf{a}_n.
```

因此，矩阵乘向量的本质是：用输入向量 $`\mathbf{x}`$ 的各个分量作为系数，对矩阵各列进行线性组合。

### 3. 列空间

矩阵各列的所有线性组合组成矩阵的**列空间**：

```math
\mathrm{Col}(A)
= \mathrm{span}\{\mathbf{a}_1,\ldots,\mathbf{a}_n\}
= \{A\mathbf{x}\mid\mathbf{x}\in\mathbb{R}^n\}.
```

每个列向量都有 $`m`$ 个分量，所以

```math
\mathrm{Col}(A)\subseteq\mathbb{R}^m.
```

列空间描述的是矩阵能够产生的全部输出，不一定等于整个 $`\mathbb{R}^m`$。

### 4. 右端向量与解的存在性

在 $`A\mathbf{x}=\mathbf{b}`$ 中，$`\mathbf{b}`$ 称为 **right-hand side（右端向量，RHS）**，它是希望矩阵产生的目标输出。

由于 $`A\mathbf{x}`$ 是矩阵各列的线性组合，所以

```math
A\mathbf{x}=\mathbf{b}\text{ 有解}
\quad\Longleftrightarrow\quad
\mathbf{b}\in\mathrm{Col}(A).
```

方程数与未知数数只能提供线索，不能单独决定某个具体方程组是否有解：

- 当 $`m>n`$ 时，方程较多，某些条件可能冲突，但特定的 $`\mathbf{b}`$ 仍可能有解。
- 当 $`m<n`$ 时，也可能因为方程冲突而无解；一旦有解，则一定存在自由变量，所以不会只有唯一解。

判断 $`\mathbf{b}`$ 能否到达之后，还需要确定哪些列真正扩大了列空间、哪些列只是已有方向的重复组合。

## 三、列空间中的独立方向、基与秩

### 1. 冗余列与线性相关

考虑贯穿后文的矩阵

```math
A=
\begin{bmatrix}
1&2&3\\
2&4&6\\
2&6&8\\
2&8&10
\end{bmatrix}
= \begin{bmatrix}
\mathbf{a}_1&\mathbf{a}_2&\mathbf{a}_3
\end{bmatrix}.
```

它的第三列满足

```math
\mathbf{a}_3=\mathbf{a}_1+\mathbf{a}_2.
```

因此第三列没有提供新的方向：

```math
\mathrm{span}\{\mathbf{a}_1,\mathbf{a}_2,\mathbf{a}_3\}
= \mathrm{span}\{\mathbf{a}_1,\mathbf{a}_2\}.
```

同时

```math
\mathbf{a}_1+\mathbf{a}_2-\mathbf{a}_3=\mathbf{0}
```

给出了一组不全为零的系数，所以三列线性相关。能够由其他列线性组合得到的列是冗余列，它可以从列空间的生成集合中去掉。

删除冗余列会改变矩阵和未知数的数量，却不会改变能够产生的输出集合。若

```math
B=
\begin{bmatrix}
\mathbf{a}_1&\mathbf{a}_2
\end{bmatrix},
```

则 $`A\neq B`$，但

```math
\mathrm{Col}(A)=\mathrm{Col}(B).
```

所以“列空间相同”不等于“方程组相同”。

### 2. 从基到维数与秩

上一节已经得到

```math
\mathrm{Col}(A)
= \mathrm{span}\{\mathbf{a}_1,\mathbf{a}_2,\mathbf{a}_3\}
= \mathrm{span}\{\mathbf{a}_1,\mathbf{a}_2\}.
```

三列都能用来生成列空间，但 $`\mathbf{a}_3`$ 可以由前两列线性组合得到，因此它是冗余的。去掉它以后，$`\mathbf{a}_1,\mathbf{a}_2`$ 仍能张成整个列空间，而且二者线性无关。

一组向量同时满足下面两个条件，就称为这个空间的一组**基**（basis）：

1. **能够张成整个空间**：空间中的每个向量都能由它们线性组合得到；
2. **线性无关**：其中没有向量可以由其他向量线性组合得到。

因此，本例中

```math
\{\mathbf{a}_1,\mathbf{a}_2\}
```

是 $`\mathrm{Col}(A)`$ 的一组基。基回答的是“哪些向量构成这个空间的一套独立生成方向”，所以**基是一组向量**。

接下来数这组基中有几个向量。基中向量的个数称为空间的**维数**（dimension）。本例的基含有两个向量，所以

```math
\dim\bigl(\mathrm{Col}(A)\bigr)=2.
```

对于矩阵 $`A`$，其**秩**（rank）定义为列空间的维数，因此

```math
\mathrm{rank}(A)
= \dim\bigl(\mathrm{Col}(A)\bigr)
=2.
```

这里不能写成“基等于秩”。三者是不同类型的对象：

| 概念 | 回答的问题 | 它是什么 | 本例 |
| --- | --- | --- | --- |
| 基 | 哪些向量是独立的生成方向？ | 一组向量 | $`\{\mathbf{a}_1,\mathbf{a}_2\}`$ |
| 维数 | 这个空间有几个独立方向？ | 一个数字 | $`2`$ |
| 秩 | 矩阵的列空间有几个独立方向？ | 一个数字 | $`2`$ |

>基属于所研究的向量空间；秩是针对矩阵而言的数字。基是一组向量；秩是这组基里有几个向量。

研究矩阵 $`A`$ 的列空间时，二者通过“维数”联系起来：

```math
\boxed{
\text{一组基所含向量的数量}
= \dim\bigl(\mathrm{Col}(A)\bigr)
= \mathrm{rank}(A)
}
```

同一个空间可以有多组不同的基。例如，下面两组向量都能作为 $`\mathbb{R}^2`$ 的基：

```math
\left\{
\begin{bmatrix}1\\0\end{bmatrix},
\begin{bmatrix}0\\1\end{bmatrix}
\right\},
\qquad
\left\{
\begin{bmatrix}1\\1\end{bmatrix},
\begin{bmatrix}1\\-1\end{bmatrix}
\right\}.
```

它们所含的具体向量不同，但每组都有两个向量。一个空间的基可以不唯一，任意一组基所含向量的数量却一定相同，因此空间的维数是确定的。

在简单矩阵中，可以直接观察哪些列是冗余的；对于一般矩阵，需要通过高斯消元找到主元列的位置。主元列的数量给出秩，而原矩阵中对应位置的列构成列空间的一组基：

```math
\boxed{
\text{主元列数量}
= \text{基向量数量}
= \dim\bigl(\mathrm{Col}(A)\bigr)
= \mathrm{rank}(A)
}
```

下一节将沿着这条关系，用高斯消元系统地识别主元列、基与秩。

## 四、用高斯消元识别主元结构

### 1. 消元的目的

高斯消元不只是机械地制造零，而是在把矩阵中隐藏的独立与依赖关系显露出来。

对矩阵做初等行变换，相当于左乘某个可逆矩阵 $`E`$：

```math
A\longmapsto EA.
```

若消元后发现

```math
E\mathbf{a}_2=cE\mathbf{a}_1,
```

则由 $`E`$ 可逆可得

```math
\mathbf{a}_2=c\mathbf{a}_1.
```

因此，行操作不会改变列之间的线性依赖关系，也不会改变齐次方程的解。

### 2. 行阶梯形与简化行阶梯形

为了使用同一个完整例子，令

```math
C=A^T
= \begin{bmatrix}
1&2&2&2\\
2&4&6&8\\
3&6&8&10
\end{bmatrix}.
```

向下消元得到行阶梯形（row echelon form）：

```math
C\longrightarrow
U=
\begin{bmatrix}
1&2&2&2\\
0&0&2&4\\
0&0&0&0
\end{bmatrix}.
```

$`U`$ 的主要特征是主元呈楼梯状排列，并且每个主元下面都是零。这里的主元位于第 $`1`$、$`3`$ 列，第二个主元仍为 $`2`$。

继续消元：

```math
R_2\leftarrow\frac{1}{2}R_2,
\qquad
R_1\leftarrow R_1-2R_2,
```

得到

```math
R=
\mathrm{rref}(C)
= \begin{bmatrix}
1&2&0&-2\\
0&0&1&2\\
0&0&0&0
\end{bmatrix}.
```

简化行阶梯形（reduced row echelon form，RREF）必须满足：

1. 所有非零行位于零行上方。
2. 每一行的主元都比上一行的主元更靠右。
3. 每个主元等于 $`1`$。
4. 每个主元所在列除主元外，其余元素全部为 $`0`$。

所以，把主元变成 $`1`$ 只是其中一步；还必须把主元上方的元素也消成零。RREF 与矩阵是否为方阵无关，任意尺寸的矩阵都可以化为 RREF，并且每个矩阵的 RREF 是唯一的。

### 3. 主元列与自由列

在 $`U`$ 或 $`R`$ 中，第 $`1`$、$`3`$ 列包含主元，因此它们是 **pivot columns（主元列）**；第 $`2`$、$`4`$ 列不包含主元，因此它们是 **free columns（自由列或非主元列）**。

主元列代表扫描列时出现了新的独立方向；非主元列可以由此前的主元列线性组合得到。这里需要区分两个现象：

| 消元结果 | 表示的含义 |
| --- | --- |
| 出现零行 | 某个方程是冗余约束 |
| 某列没有主元 | 该列没有提供新的独立方向 |

某个候选主元位置变成零，并不表示“这一行没有意义”；它表示在该列中找不到新的主元。

### 4. 用消元结果寻找列空间的基

主元列的编号由消元后的矩阵确定，但 $`\mathrm{Col}(C)`$ 的基必须取自**原矩阵 $`C`$ 的对应列**。

这是因为行操作保持主元位置、秩和列之间的依赖关系，却通常会改变列向量本身以及它们实际张成的列空间：

```math
\mathrm{Col}(EC)=E\mathrm{Col}(C),
```

一般并不等于 $`\mathrm{Col}(C)`$。

本例的主元列编号是 $`1,3`$，所以 $`\mathrm{Col}(C)`$ 的一组基是

```math
\left\{
\begin{bmatrix}1\\2\\3\end{bmatrix},
\begin{bmatrix}2\\6\\8\end{bmatrix}
\right\},
```

而不是直接取 $`R`$ 的第 $`1`$、$`3`$ 列作为原列空间的基。

主元的数量就是秩：

```math
r=主元个数=rank(C)=2.
```

主元位置不仅筛选出列空间的独立方向，还会把未知数分成主元变量和自由变量。自由变量正是进入零空间的入口。

## 五、从自由变量到零空间

### 1. 齐次方程与零空间

当右端向量为零时，方程变成

```math
C\mathbf{x}=\mathbf{0}.
```

这叫作**齐次线性方程组**（homogeneous system）。无论 $`C`$ 是什么矩阵，$`\mathbf{x}=\mathbf{0}`$ 都是一个解，称为**平凡解**（trivial solution）。矩阵 $`C`$ 的结构决定的是除平凡解外是否还有非零解。

所有满足齐次方程的输入组成 $`C`$ 的零空间：

```math
N(C)=\{\mathbf{x}\in\mathbb{R}^n\mid C\mathbf{x}=\mathbf{0}\}.
```

从列向量角度看，零空间记录的是所有能让矩阵各列线性组合成零向量的系数组合。接下来利用消元后的主元结构求出这个空间。

### 2. 从列的位置对应到变量

在方程

```math
C\mathbf{x}=\mathbf{0}
```

中，每一列对应一个未知数：

| 列 | 第1列 | 第2列 | 第3列 | 第4列 |
| --- | --- | --- | --- | --- |
| 类型 | 主元列 | 自由列 | 主元列 | 自由列 |
| 变量 | $`x_1`$ | $`x_2`$ | $`x_3`$ | $`x_4`$ |
| 变量类型 | 主元变量 | 自由变量 | 主元变量 | 自由变量 |

自由变量（F free）可以独立指定；主元变量（I identity）不是固定常数，而是由自由变量决定。

### 3. 行操作与右端向量

研究一般方程 $`C\mathbf{x}=\mathbf{b}`$ 时，必须对增广矩阵

```math
[C\mid\mathbf{b}]
```

整体执行相同的行操作。

研究零空间时，方程按定义是

```math
C\mathbf{x}=\mathbf{0}.
```

零向量经过任何行操作仍然是零向量，所以右侧的零列通常省略，只对系数矩阵写出

```math
C\longrightarrow U\longrightarrow R.
```

这些行操作可逆，因此

```math
N(C)=N(U)=N(R).
```

这里的 $`R`$ 只是系数矩阵的 RREF，并不是增广矩阵。

### 4. 从 RREF 读出齐次方程的解

由

```math
R=
\begin{bmatrix}
1&2&0&-2\\
0&0&1&2\\
0&0&0&0
\end{bmatrix}
```

可得

```math
\begin{cases}
x_1+2x_2-2x_4=0,\\
x_3+2x_4=0.
\end{cases}
```

令自由变量

```math
x_2=s,
\qquad
x_4=t,
```

则主元变量为

```math
x_1=-2s+2t,
\qquad
x_3=-2t.
```

因此所有解是

```math
\mathbf{x}
= \begin{bmatrix}
-2s+2t\\
s\\
-2t\\
t
\end{bmatrix}
= s
\begin{bmatrix}
-2\\1\\0\\0
\end{bmatrix}
+
t
\begin{bmatrix}
2\\0\\-2\\1
\end{bmatrix}.
```

### 5. 特殊解与零空间的基

零空间定义为

```math
N(C)=\{\mathbf{x}\in\mathbb{R}^4\mid C\mathbf{x}=\mathbf{0}\}.
```

分别令一个自由变量为 $`1`$、其余自由变量为 $`0`$，得到 **special solutions（特殊解）**：

```math
\mathbf{v}_1=
\begin{bmatrix}-2\\1\\0\\0\end{bmatrix},
\qquad
\mathbf{v}_2=
\begin{bmatrix}2\\0\\-2\\1\end{bmatrix}.
```

所有解都是这两个特殊解的线性组合，所以

```math
N(C)
= \mathrm{span}
\left\{
\mathbf{v}_1,\mathbf{v}_2
\right\}.
```

$`\mathbf{v}_1,\mathbf{v}_2`$ 构成零空间的一组基。参数 $`s,t`$ 只是任意实数系数；自由变量有几个，就会得到几个基本的特殊解方向。

零空间中的向量记录的是矩阵各列之间的依赖关系。例如

```math
C\mathbf{v}_1=\mathbf{0}
\quad\Longleftrightarrow\quad
-2\mathbf{c}_1+\mathbf{c}_2=\mathbf{0},
```

而

```math
C\mathbf{v}_2=\mathbf{0}
\quad\Longleftrightarrow\quad
2\mathbf{c}_1-2\mathbf{c}_3+\mathbf{c}_4=\mathbf{0}.
```

因此，零空间可以理解为所有“让矩阵各列抵消成零的系数组合”。

### 6. RREF 的分块形式

在 RREF 中，把主元列与自由列分别集中排列，并忽略底部零行，可以把非零部分写成

```math
R=\begin{bmatrix}I&F\end{bmatrix}.
```

$`I`$ 是主元列形成的单位矩阵块，$`F`$ 是自由列组成的矩阵块。真正可以自由选择的是与 $`F`$ 对应的变量 $`\mathbf{x}_F`$，不是矩阵 $`F`$ 本身。

把变量也按主元变量、自由变量分组：

```math
\mathbf{x}
= \begin{bmatrix}
\mathbf{x}_P\\
\mathbf{x}_F
\end{bmatrix}.
```

齐次方程成为

```math
\begin{bmatrix}I&F\end{bmatrix}
\begin{bmatrix}
\mathbf{x}_P\\
\mathbf{x}_F
\end{bmatrix}
= \mathbf{0}.
```

于是

```math
\mathbf{x}_P+F\mathbf{x}_F=\mathbf{0},
```

所以

```math
\mathbf{x}_P=-F\mathbf{x}_F.
```

负号只来自移项。整个解可以写为

```math
\begin{bmatrix}
\mathbf{x}_P\\
\mathbf{x}_F
\end{bmatrix}
= \begin{bmatrix}
-F\\
I
\end{bmatrix}
\mathbf{x}_F.
```

矩阵

```math
\begin{bmatrix}-F\\I\end{bmatrix}
```

的各列就是零空间的特殊解，也构成零空间的一组基。

这个紧凑公式有一个重要前提：列和变量已经按“主元在前、自由变量在后”重新排列。当前例子的顺序是

```math
x_1,\ x_2,\ x_3,\ x_4
\quad\longrightarrow\quad
x_1,\ x_3\mid x_2,\ x_4.
```

如果主元列和自由列在原矩阵中交错排列，最后必须把解向量的分量恢复到原变量顺序，不能直接机械地上下拼接。

自由变量的数量就是零空间的维数。这个数量与主元数之间的关系，把零空间、秩以及转置后的矩阵联系起来。

## 六、秩、零度与转置

### 1. 秩与自由变量数量

若

```math
C\in\mathbb{R}^{m\times n},
\qquad
\mathrm{rank}(C)=r,
```

则主元变量有 $`r`$ 个，自由变量有 $`n-r`$ 个。因此

```math
\mathrm{nullity}(C)
= \dim N(C)
= n-r.
```

这就是秩－零度定理：

```math
\mathrm{rank}(C)+\mathrm{nullity}(C)=n.
```

当前例子中 $`n=4,r=2`$，所以

```math
\dim N(C)=4-2=2.
```

没有自由变量时，齐次系统只有平凡解；只要存在自由变量，齐次系统就有无穷多个解。

### 2. 转置与秩

矩阵转置会交换行和列，但不会改变秩：

```math
\mathrm{rank}(A)=\mathrm{rank}(A^T).
```

其本质是行秩等于列秩。对前面的矩阵，

```math
\mathrm{rank}(A)
= \mathrm{rank}(A^T)
=2.
```

不过，转置前后的零空间通常不同。因为

```math
A:\mathbb{R}^3\to\mathbb{R}^4,
\qquad
A^T:\mathbb{R}^4\to\mathbb{R}^3,
```

所以

```math
N(A)\subseteq\mathbb{R}^3,
\qquad
N(A^T)\subseteq\mathbb{R}^4.
```

本例中

```math
\dim N(A)=3-2=1,
\qquad
\dim N(A^T)=4-2=2.
```

秩相同不代表零空间相同。

到这里，列空间已经回答目标输出是否能够到达，零空间则描述不会改变输出的输入变化。二者合在一起，就能写出非齐次方程的全部解。

## 七、用列空间与零空间描述完整解

### 1. 特解与零空间解

假设 $`A\mathbf{x}=\mathbf{b}`$ 有解，并且已经找到其中一个解 $`\mathbf{x}_p`$：

```math
A\mathbf{x}_p=\mathbf{b}.
```

$`\mathbf{x}_p`$ 称为 **particular solution（特解）**。再取零空间中的任意向量 $`\mathbf{x}_n`$：

```math
\mathbf{x}_n\in N(A),
\qquad
A\mathbf{x}_n=\mathbf{0}.
```

利用线性性可得

```math
A(\mathbf{x}_p+\mathbf{x}_n)
=A\mathbf{x}_p+A\mathbf{x}_n
=\mathbf{b}+\mathbf{0}
=\mathbf{b}.
```

因此，在特解上加任意零空间向量，仍然是原方程的解。这里能加的不是任意向量，而只能是 $`N(A)`$ 中的向量。

### 2. 为什么这包含全部解

若 $`\mathbf{x}`$ 是 $`A\mathbf{x}=\mathbf{b}`$ 的任意一个解，那么

```math
A\mathbf{x}=\mathbf{b},
\qquad
A\mathbf{x}_p=\mathbf{b}.
```

两式相减：

```math
A(\mathbf{x}-\mathbf{x}_p)=\mathbf{0}.
```

所以

```math
\mathbf{x}-\mathbf{x}_p\in N(A).
```

令 $`\mathbf{x}_n=\mathbf{x}-\mathbf{x}_p`$，便得到

```math
\boxed{
\mathbf{x}=\mathbf{x}_p+\mathbf{x}_n,
\qquad
\mathbf{x}_n\in N(A)
}.
```

这说明所有解，而且只有所有解，都属于集合

```math
\boxed{
\{\mathbf{x}:A\mathbf{x}=\mathbf{b}\}
=\mathbf{x}_p+N(A)
}.
```

因此可以概括为

```math
\text{complete solution（通解）}
= \text{particular solution（特解）}
+
\text{null-space solution（零空间解）}.
```

### 3. 自由变量在完整解中的作用

求一个方便的特解时，通常先把所有自由变量设为 $`0`$，再计算主元变量。这样得到的 $`\mathbf{x}_p`$ 只负责产生目标输出：

```math
A\mathbf{x}_p=\mathbf{b}.
```

自由变量能够产生的全部变化集中在 $`\mathbf{x}_n`$ 中，而这些变化不会影响输出：

```math
A\mathbf{x}_n=\mathbf{0}.
```

从列向量视角看，$`\mathbf{x}_p`$ 是一套把矩阵各列组合成 $`\mathbf{b}`$ 的系数；$`\mathbf{x}_n`$ 是一套让各列互相抵消成零的系数。两套系数相加后，输出仍然是 $`\mathbf{b}`$。

### 4. 解集的几何形状

零空间 $`N(A)`$ 一定经过原点，是一个子空间。若方程有一个特解 $`\mathbf{x}_p`$，则非齐次方程的解集是

```math
\mathbf{x}_p+N(A),
```

也就是把零空间整体平移到经过 $`\mathbf{x}_p`$ 的位置。

- 若 $`N(A)`$ 是一条直线，解集是经过 $`\mathbf{x}_p`$、与该直线平行的直线。
- 若 $`N(A)`$ 是一个平面，解集是经过 $`\mathbf{x}_p`$、与该平面平行的平面。

这种平移后的集合称为**仿射空间**（affine space）。当 $`\mathbf{b}\neq\mathbf{0}`$ 时，解集不包含原点，因此不是子空间；当 $`\mathbf{b}=\mathbf{0}`$ 时，解集就是 $`N(A)`$，它是子空间。

完整解公式把问题分成了两步：列空间决定特解是否存在，零空间决定能否在特解之外继续变化。接下来只需比较 $`r,m,n`$，就能统一判断解的数量。

## 八、用秩判断解的数量

### 1. 两个核心判断

设

```math
A\in\mathbb{R}^{m\times n},
\qquad
r=\mathrm{rank}(A).
```

其中 $`m`$ 是行数和方程数，$`n`$ 是列数和未知数数，$`r`$ 是主元数量。

解的结构由两个条件分别控制：

```math
r=m
\quad\Longleftrightarrow\quad
\text{每一行都有主元}
\quad\Longleftrightarrow\quad
\mathrm{Col}(A)=\mathbb{R}^m.
```

所以 $`r=m`$ 保证 $`A\mathbf{x}=\mathbf{b}`$ 对每个 $`\mathbf{b}\in\mathbb{R}^m`$ 都有解，控制的是**存在性**。

另一方面，

```math
r=n
\quad\Longleftrightarrow\quad
\text{每一列都有主元}
\quad\Longleftrightarrow\quad
N(A)=\{\mathbf{0}\}.
```

所以 $`r=n`$ 保证没有自由变量，方程只要有解就必然唯一，控制的是**唯一性**。

### 2. 四种典型情形

| 秩与尺寸的关系 | 分组后的 RREF 结构 | 自由变量 | 是否对每个 $`\mathbf{b}`$ 有解 | 某个 $`A\mathbf{x}=\mathbf{b}`$ 的解数 |
| --- | --- | ---: | --- | --- |
| $`r=m=n`$ | $`I`$ | $`0`$ | 是 | 恰好一个 |
| $`r=n<m`$ | $`\begin{bmatrix}I\\0\end{bmatrix}`$ | $`0`$ | 不一定 | $`0`$ 或 $`1`$ |
| $`r=m<n`$ | $`\begin{bmatrix}I&F\end{bmatrix}`$ | $`n-r>0`$ | 是 | 无穷多个 |
| $`r<m`$ 且 $`r<n`$ | $`\begin{bmatrix}I&F\\0&0\end{bmatrix}`$ | $`n-r>0`$ | 不一定 | $`0`$ 或无穷多个 |

表中的 $`I`$ 与 $`F`$ 是把主元列和自由列重新分组后的抽象结构；它们在原矩阵中可能交错排列。

#### 满秩方阵：$`r=m=n`$

每一行和每一列都有主元，因此所有 $`\mathbf{b}`$ 都可到达，同时没有自由变量。于是每个 $`\mathbf{b}`$ 都有唯一解。

#### 满列秩高矩阵：$`r=n<m`$

每一列都有主元，所以没有自由变量；但不是每一行都有主元，因此列空间不能填满 $`\mathbb{R}^m`$。某个方程可能无解，也可能有唯一解，但不可能有无穷多个解。

#### 满行秩宽矩阵：$`r=m<n`$

每一行都有主元，所以每个 $`\mathbf{b}`$ 都有解；同时 $`n-r>0`$，存在自由变量，因此每个 $`\mathbf{b}`$ 都有无穷多个解。

#### 行列都不满秩：$`r<m`$ 且 $`r<n`$

不是每个 $`\mathbf{b}`$ 都在列空间中，所以某些方程无解；一旦有解，又会因为存在自由变量而有无穷多个解。

### 3. 满秩方阵与唯一解

当

```math
r=m=n
```

时，$`A`$ 是满秩方阵。这里“没有任何”指的是没有自由变量、没有非零的零空间解，而不是没有特解。

因为 $`r=n`$，所以

```math
N(A)=\{\mathbf{0}\}.
```

完整解公式退化为

```math
\mathbf{x}=\mathbf{x}_p+\mathbf{0}=\mathbf{x}_p.
```

同时，因为 $`r=m`$，任意 $`\mathbf{b}`$ 都在列空间中，所以特解总是存在。因此每个右端向量都有且只有一个解。

满秩方阵的 RREF 是单位矩阵：

```math
\mathrm{rref}(A)=I.
```

注意，这是说 $`A`$ 的 RREF 等于 $`I`$，并不是说原矩阵 $`A`$ 本身一定等于 $`I`$。

### 4. 可逆矩阵

对 $`n\times n`$ 方阵，满秩还意味着 $`A`$ 可逆。此时 $`A^{-1}`$ 存在，而且

```math
A\mathbf{x}=\mathbf{b}
\quad\Longrightarrow\quad
\mathbf{x}=A^{-1}\mathbf{b}.
```

例如

```math
A=
\begin{bmatrix}
1&2\\
3&1
\end{bmatrix}
```

满足

```math
\det(A)=1\cdot1-2\cdot3=-5\neq0,
```

所以它可逆，并且对任意 $`\mathbf{b}\in\mathbb{R}^2`$ 都有唯一解。

对于 $`n\times n`$ 方阵，下面这些说法彼此等价：

```math
\mathrm{rank}(A)=n
\Longleftrightarrow
\text{每行、每列都有主元}
\Longleftrightarrow
\mathrm{rref}(A)=I
```

```math
\Longleftrightarrow
N(A)=\{\mathbf{0}\}
\Longleftrightarrow
\mathrm{Col}(A)=\mathbb{R}^n
\Longleftrightarrow
A^{-1}\text{ 存在}
```

```math
\Longleftrightarrow
A\mathbf{x}=\mathbf{b}
\text{ 对每个 }\mathbf{b}\in\mathbb{R}^n
\text{ 都有唯一解}.
```

### 5. 解的数量为什么只有三种

线性方程组的解数只能是

```math
0,\qquad1,\qquad\infty.
```

如果 $`\mathbf{x}^{(1)}`$ 和 $`\mathbf{x}^{(2)}`$ 是两个不同的解，那么

```math
A(\mathbf{x}^{(1)}-\mathbf{x}^{(2)})=\mathbf{0}.
```

令

```math
\mathbf{v}=\mathbf{x}^{(1)}-\mathbf{x}^{(2)}\neq\mathbf{0},
```

则对任意 $`c\in\mathbb{R}`$，

```math
A(\mathbf{x}^{(1)}+c\mathbf{v})
=\mathbf{b}+c\mathbf{0}
=\mathbf{b}.
```

因此，只要出现两个不同的解，就会沿着非零零空间方向产生无穷多个解，不可能恰好只有 $`2`$ 个、$`3`$ 个或其他有限多个解。

## 九、当前课程的统一逻辑

### 1. 输出端与输入端

对于 $`C\in\mathbb{R}^{m\times n}`$：

```math
\mathrm{Col}(C)\subseteq\mathbb{R}^m,
\qquad
N(C)\subseteq\mathbb{R}^n.
```

二者从不同方向描述同一个矩阵：

| 空间 | 核心问题 |
| --- | --- |
| 列空间 $`\mathrm{Col}(C)`$ | 哪些输出 $`\mathbf{b}`$ 能由 $`C`$ 产生？ |
| 零空间 $`N(C)`$ | 哪些输入 $`\mathbf{x}`$ 会被 $`C`$ 送到零向量？ |

因此：

```math
C\mathbf{x}=\mathbf{b}\text{ 是否有解}
\quad\text{由列空间决定；}
```

```math
C\mathbf{x}=\mathbf{b}\text{ 若有解，是否唯一}
\quad\text{由零空间决定。}
```

若 $`\mathbf{x}_p`$ 是 $`C\mathbf{x}=\mathbf{b}`$ 的一个特解，那么所有解都可以写成

```math
\mathbf{x}=\mathbf{x}_p+\mathbf{x}_n,
\qquad
\mathbf{x}_n\in N(C).
```

### 2. 整堂课的逻辑链

```math
\text{线性组合}
\longrightarrow
\text{张成子空间}
\longrightarrow
\text{矩阵各列张成列空间}
```

```math
A\mathbf{x}=\mathbf{b}\text{ 有解}
\Longleftrightarrow
\mathbf{b}\in\mathrm{Col}(A)
```

```math
\text{高斯消元}
\longrightarrow
\text{主元位置}
\longrightarrow
\text{秩与独立列}
\longrightarrow
\text{列空间的基}
```

```math
\text{非主元列}
\longrightarrow
\text{自由变量}
\longrightarrow
\text{特殊解}
\longrightarrow
\text{零空间的基}
```

```math
\mathbf{b}\in\mathrm{Col}(A)
\longrightarrow
\text{找到一个特解 }\mathbf{x}_p
\longrightarrow
\mathbf{x}=\mathbf{x}_p+\mathbf{x}_n,
\quad \mathbf{x}_n\in N(A)
```

```math
r=m\Longrightarrow\text{每个右端向量都有解},
\qquad
r=n\Longrightarrow\text{有解时解唯一}.
```

最终，列空间描述矩阵能到达哪里，零空间描述矩阵会丢失哪些输入方向；特解负责到达目标输出，零空间负责描述不改变输出的全部自由变化，而主元和秩把这些结构联系在一起。

### 3. 术语对照

| 英文 | 中文 | 本节含义 |
| --- | --- | --- |
| scalar | 标量 | 与向量相乘的普通数 |
| span | 张成 | 一组向量所有线性组合的集合 |
| column space | 列空间 | 矩阵所有可能输出的集合 |
| right-hand side | 右端向量 | $`A\mathbf{x}=\mathbf{b}`$ 中的目标输出 $`\mathbf{b}`$ |
| pivot | 主元 | 阶梯形矩阵中每个非零行的领先非零元素 |
| pivot column | 主元列 | 包含主元、提供独立方向的列 |
| free column | 自由列或非主元列 | 不包含主元的列 |
| pivot variable | 主元变量 | 由方程和自由变量决定的变量 |
| free variable | 自由变量 | 可以独立选取的参数 |
| row echelon form | 行阶梯形 | 主元下方为零的阶梯形矩阵 |
| RREF | 简化行阶梯形 | 主元为 $`1`$，且主元列其余位置全为 $`0`$ |
| rank | 秩 | 主元数量，也是列空间的维数 |
| null space | 零空间 | 所有满足 $`A\mathbf{x}=\mathbf{0}`$ 的输入组成的空间 |
| nullity | 零度 | 零空间的维数，即自由变量数量 |
| special solution | 特殊解 | 依次把一个自由变量设为 $`1`$、其余设为 $`0`$ 得到的解 |
| particular solution | 特解 | $`A\mathbf{x}=\mathbf{b}`$ 的某一个具体解 $`\mathbf{x}_p`$ |
| complete solution | 通解 | $`\mathbf{x}_p+N(A)`$，即非齐次方程的全部解 |
| affine space | 仿射空间 | 子空间经过平移后得到的集合 |
| full rank | 满秩 | 秩达到 $`\min(m,n)`$ |
| invertible matrix | 可逆矩阵 | 存在逆矩阵 $`A^{-1}`$ 的方阵 |
