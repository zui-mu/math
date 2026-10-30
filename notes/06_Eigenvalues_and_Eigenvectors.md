# Eigenvalues and Eigenvectors

## 本章主线与范围

上一章用 determinant 判断方阵是否可逆。本章把这个工具用于一个新问题：**给定一个线性变换，能否找到一些非零向量，使矩阵在这些向量上只表现为数乘？**

```math
\boxed{Av=\lambda v,\qquad v\neq\mathbf0}.
```

本章以实矩阵和实向量的几何直觉为主；讨论特征多项式的全部根时，允许在复数域中求解。实矩阵不一定有实特征值，不能默认每个矩阵在实空间中都有不变直线。

```text
已知方阵 A，寻找未知的数 λ 和非零向量 v
    -> Av=λv：A 在某个方向上只做数乘
    -> λv=(λI)v，得到 (A-λI)v=0
    -> 非零解存在 <=> A-λI singular <=> det(A-λI)=0
    -> 先求 λ，再求 N(A-λI)，得到特征空间
    -> 用具体矩阵区分：变换规则、输入向量、列空间与零空间
    -> 投影矩阵：固定的方向对应 1，被消去的方向对应 0
    -> 特征多项式：特征值之和是 trace，乘积是 determinant
    -> 区分重复特征根与独立特征方向，判断能否对角化
    -> 用特征向量基表示初始向量，解 u_(k+1)=A u_k
    -> 斐波那契：递推关系 -> 二维状态 -> 特征分解 -> 通式
```

相关前置笔记：[基本子空间](01_向量空间与四个基本子空间.md)、[正交投影](04_Projections_onto_Subspaces.md)、[行列式](05_Determinants.md)。

> **本章不是把零空间换个名字：我们先寻找哪些 $`\lambda`$ 能让 $`A-\lambda I`$ 有非零零空间，再研究这些非零解怎样描述原矩阵 $`A`$ 的作用。**

## 一、这个公式在问什么：已知、未知与几何意义

### 1. 已知矩阵，寻找特征对

通常的求解任务是：

| 对象 | 已知还是未知 | 含义与尺寸 |
| --- | --- | --- |
| $`A`$ | 已知 | $`n\times n`$ 方阵，给定的线性变换规则 |
| $`\lambda`$ | 未知 | eigenvalue（特征值），一个标量 |
| $`v`$ | 未知 | eigenvector（特征向量），$`n\times1`$ 非零列向量 |
| $`Av`$ | 由 $`A,v`$ 决定 | 变换后的向量，仍有 $`n`$ 个分量 |

**先有矩阵 $`A`$，再寻找满足关系的 $`(\lambda,v)`$，不是先随便选一个 $`\lambda`$，就一定能找到对应的非零向量。**

与解 $`Ax=b`$ 相比，区别是：这里没有单独给定目标向量 $`b`$；输出必须恰好等于同一个输入向量 $`v`$ 的某个倍数。

采用列向量约定时，$`Av`$ 是**矩阵乘向量**，不是把向量放在左边乘矩阵。只有输入、输出处于同一个空间，才可以直接比较 $`Av`$ 与 $`v`$，所以通常讨论方阵的特征值。

### 2. “保持方向”应理解为保持同一条直线

当 $`\lambda`$ 为实数时：

| 特征值 | 对特征向量的作用 |
| --- | --- |
| $`\lambda>1`$ | 同向，长度放大 $`\lambda`$ 倍 |
| $`0<\lambda<1`$ | 同向，长度缩小 |
| $`\lambda=1`$ | 向量保持不变 |
| $`\lambda<0`$ | 反向，长度乘 $`\lvert\lambda\rvert`$ |
| $`\lambda=0`$ | 非零输入被送到零向量 |

所以更准确的几何描述是

```math
A\bigl(\mathrm{span}\{v\}\bigr)\subseteq\mathrm{span}\{v\}.
```

这条过原点的直线是一个 **invariant subspace（不变子空间）**。负特征值允许方向反转；零特征值允许输出为零，不能说零输出仍有一个几何方向。

若把任意实 $`\lambda`$ 的统一缩放说成“只拉长”，也不准确：它还可能缩短、反向或归零。

### 3. 为什么排除零向量，为什么解不唯一

因为对每个标量 $`\lambda`$ 都有

```math
A\mathbf0=\lambda\mathbf0=\mathbf0.
```

如果允许 $`v=\mathbf0`$，这个等式对 $`\lambda`$ 就没有筛选作用。因此必须要求 $`v\neq\mathbf0`$。

若 $`Av=\lambda v`$ 且 $`c\neq0`$，则

```math
A(cv)=cAv=\lambda(cv).
```

所以特征向量的**非零倍数**仍是同一特征值的特征向量。通常求的是整个特征空间的一组基，而不是唯一的一根有固定长度的向量。

> **零向量不是特征向量；特征空间却必须包含零向量，因为特征空间是一个线性子空间。**

## 二、从原方程到齐次方程：单位矩阵为何出现

### 1. 数乘、矩阵缩放和行列式性质是三件事

例如 $`v=(x,y)^T`$，则

```math
\lambda v=
\begin{bmatrix}\lambda x\\\lambda y\end{bmatrix}.
```

标量要乘向量的每个分量。为了把它写成矩阵乘法，使用

```math
I_2v=v,\qquad
\lambda I_2=
\begin{bmatrix}\lambda&0\\0&\lambda\end{bmatrix}.
```

因此

```math
\lambda v=(\lambda I_n)v.
```

不是为了“增加一个新变换”，而是利用单位矩阵的作用不变，把两边改写成可相减的同尺寸矩阵作用。

与上一章的规则对照：

| 运算 | 正确含义 |
| --- | --- |
| $`\lambda v`$ | 向量每个分量乘 $`\lambda`$ |
| $`\lambda A`$ | 矩阵每个元素乘 $`\lambda`$ |
| $`\det(\lambda A)`$ | 对 $`n\times n`$ 方阵，等于 $`\lambda^n\det(A)`$ |
| 某一行单独乘 $`\lambda`$ | 得到一个行缩放后的矩阵，其 determinant 乘 $`\lambda`$ |

行列式中“从某一行提出一个因子”的规则，不能用来把向量数乘改成只缩放一个分量，也不能据此声称缩放一行与缩放整个矩阵相同。

### 2. 只从主对角线减去特征值

从定义出发：

```math
\begin{aligned}
Av&=\lambda v,\\
Av-(\lambda I_n)v&=\mathbf0,\\
(A-\lambda I_n)v&=\mathbf0.
\end{aligned}
```

若

```math
A=\begin{bmatrix}a&b\\c&d\end{bmatrix},
```

则

```math
\boxed{
A-\lambda I_2=
\begin{bmatrix}
a-\lambda&b\\
c&d-\lambda
\end{bmatrix}
}.
```

只改主对角线，非对角线不变。$`(A-I)v=0`$ 是其中 $`\lambda=1`$ 的特殊情况，不是一般的特征向量方程。

### 3. 为什么不能把两边的 $`v`$ 约掉

从

```math
Av=(\lambda I_n)v
```

不能推出 $`A=\lambda I_n`$。它只表示两个矩阵对**这个向量**的作用相同，不表示对所有输入都相同。

具体反例：

```math
A=\begin{bmatrix}2&0\\0&3\end{bmatrix},
\qquad
v=\begin{bmatrix}1\\0\end{bmatrix}.
```

虽然

```math
Av=\begin{bmatrix}2\\0\end{bmatrix}=2v=(2I_2)v,
```

但 $`A\neq2I_2`$，因为它对另一个方向 $`(0,1)^T`$ 的作用是乘 3。

**只有两个线性变换对一组基的每个向量都作用相同，才能据此断定矩阵相同。** 单个非零列向量不能像非零实数那样随意约去。

## 三、求解逻辑：特征方程筛选数，零空间给出向量

### 1. 为什么 determinant 必须等于零

把 $`M=A-\lambda I_n`$ 暂时看成一个整体。我们需要的不是任意解，而是

```math
Mv=\mathbf0,\qquad v\neq\mathbf0.
```

上一章已经证明：方阵的齐次方程有非零解，当且仅当矩阵不可逆。因此

```math
\boxed{
\begin{aligned}
\lambda\text{ 是 }A\text{ 的特征值}
&\Longleftrightarrow N(A-\lambda I_n)\neq\{\mathbf0\}\\
&\Longleftrightarrow \mathrm{rank}(A-\lambda I_n)<n\\
&\Longleftrightarrow A-\lambda I_n\text{ singular}\\
&\Longleftrightarrow \det(A-\lambda I_n)=0.
\end{aligned}
}
```

这里“行列式为零”是从定义推导出的**必要且充分条件**，不是人为把一个本来不为零的行列式强行设成零。

必须区分：

- 原矩阵 $`A`$ 不一定 singular；
- 对某个特征值 $`\lambda`$ 而言，移位后的矩阵 $`A-\lambda I_n`$ 一定 singular。

例如可逆矩阵也有特征值，只是不能有零特征值。

> **不要漏掉 $`\det`$：求的是 $`\det(A-\lambda I_n)=0`$，不是 $`A-\lambda I_n=0`$。前者只要求移位矩阵不可逆；后者要求它的每个元素都为零，相当于 $`A=\lambda I_n`$，是强得多的特殊情况。**

### 2. 特征空间就是移位矩阵的零空间

对任意标量 $`\lambda`$，记

```math
E_\lambda=N(A-\lambda I_n)
=\{v:Av=\lambda v\}.
```

当 $`\lambda`$ 是特征值时，$`E_\lambda`$ 称为对应的 **eigenspace（特征空间）**。其中的每一个非零向量都是特征向量。

若 $`\lambda`$ 不是特征值，则这个集合只有 $`\mathbf0`$，没有特征向量。若已经知道 $`\lambda`$ 是特征值，则

```math
\boxed{
\dim E_\lambda=n-\mathrm{rank}(A-\lambda I_n)\geq1
}.
```

特征空间不一定是一条直线。例如 $`A=2I_2`$ 时，

```math
E_2=N(2I_2-2I_2)=N(0)=\mathbb R^2.
```

此时所有非零二维向量都是特征值 2 对应的特征向量，不能只用一根向量的倍数描述全部解。

### 3. 可直接执行的求解流程

1. 写出 $`A-\lambda I_n`$，只在主对角线减去 $`\lambda`$。
2. 解 **characteristic equation（特征方程）** $`\det(A-\lambda I_n)=0`$，得到特征值。
3. 对每个不同的特征值分别代回，求 $`N(A-\lambda I_n)`$ 的一组基。
4. 用自由参数写出整个特征空间；取参数不全为零，得到全部特征向量。
5. 把基向量代回 $`Av=\lambda v`$，检查计算结果。

第二步解的是关于标量 $`\lambda`$ 的多项式方程；第三步解的是系数已经确定的齐次线性方程组。**先筛选哪些数允许非零解，再求具体的非零解。**

求 $`A-\lambda I`$ 的零空间可以正常做行消元，但一般不能先把 $`A`$ 消元成 RREF，再用 RREF 的特征值代替原矩阵的特征值：一般行操作不保持特征值。

## 四、完整算例：可逆矩阵也能产生 singular 的移位矩阵

设

```math
A=\begin{bmatrix}2&1\\1&2\end{bmatrix}.
```

它本身 $`\det(A)=3\neq0`$。下面仍能找到特征值，因为需要 singular 的是 $`A-\lambda I`$，不是 $`A`$。

### 1. 先求特征值

```math
A-\lambda I_2=
\begin{bmatrix}
2-\lambda&1\\
1&2-\lambda
\end{bmatrix}.
```

```math
\begin{aligned}
\det(A-\lambda I_2)
&=(2-\lambda)^2-1\\
&=\lambda^2-4\lambda+3\\
&=(\lambda-3)(\lambda-1).
\end{aligned}
```

所以

```math
\boxed{\lambda_1=3,\qquad\lambda_2=1}.
```

### 2. 对 $`\lambda=3`$ 求零空间

```math
A-3I_2=
\begin{bmatrix}-1&1\\1&-1\end{bmatrix}.
```

令 $`v=(x,y)^T`$，方程组为

```math
-x+y=0,\qquad x-y=0.
```

这两个方程只给出一个独立约束 $`y=x`$。令 $`x=t`$，得到

```math
E_3=
\left\{t\begin{bmatrix}1\\1\end{bmatrix}:t\in\mathbb R\right\}
=\mathrm{span}\left\{\begin{bmatrix}1\\1\end{bmatrix}\right\}.
```

对应特征向量是其中 $`t\neq0`$ 的向量。检查：

```math
A\begin{bmatrix}1\\1\end{bmatrix}
=\begin{bmatrix}3\\3\end{bmatrix}
=3\begin{bmatrix}1\\1\end{bmatrix}.
```

### 3. 对 $`\lambda=1`$ 求零空间

```math
A-I_2=\begin{bmatrix}1&1\\1&1\end{bmatrix}.
```

条件为 $`x+y=0`$，所以

```math
E_1=
\mathrm{span}\left\{\begin{bmatrix}1\\-1\end{bmatrix}\right\}.
```

检查：

```math
A\begin{bmatrix}1\\-1\end{bmatrix}
=\begin{bmatrix}1\\-1\end{bmatrix}.
```

沿 $`(1,1)^T`$ 的方向放大 3 倍，沿 $`(1,-1)^T`$ 的方向不变。注意：这里的 $`E_1=N(A-I)`$，不是 $`N(A)`$；因为 $`A`$ 可逆，$`N(A)=\{\mathbf0\}`$。

## 五、重点例子：全 1 矩阵怎样作用于二维空间

本节使用另一个矩阵，并记为 $`J`$，避免与上一节混淆：

```math
J=\begin{bmatrix}1&1\\1&1\end{bmatrix}.
```

### 1. 二维向量、矩阵的两列与两条直线

$`J`$ 把 $`\mathbb R^2`$ 映射到 $`\mathbb R^2`$。输入向量

```math
u=\begin{bmatrix}x\\y\end{bmatrix}
```

只有两个坐标：第一分量是横坐标，第二分量是纵坐标。矩阵的四个数不是一个“四维点”的坐标，而是编码了二维线性变换。

矩阵各列分别是标准基向量的输出：

```math
Je_1=J\begin{bmatrix}1\\0\end{bmatrix}
=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
Je_2=J\begin{bmatrix}0\\1\end{bmatrix}
=\begin{bmatrix}1\\1\end{bmatrix}.
```

所以对任意输入，

```math
\boxed{
Ju=xJe_1+yJe_2
=\begin{bmatrix}x+y\\x+y\end{bmatrix}
}.
```

所有输出的两个分量相等，因而落在直线 $`y=x`$ 上；而且每个 $`(s,s)^T`$ 都能由输入 $`(s,0)^T`$ 产生。所以

```math
\mathrm{Col}(J)
=\left\{\begin{bmatrix}s\\s\end{bmatrix}:s\in\mathbb R\right\}
=\mathrm{span}\left\{\begin{bmatrix}1\\1\end{bmatrix}\right\}.
```

这里必须区分：**矩阵表示变换规则，一个向量表示一根箭头；由该向量的所有倍数组成的 span 才是一整条过原点的直线。**

### 2. 零空间为什么是另一条直线

求 $`N(J)`$ 就是解

```math
Ju=\mathbf0
\Longleftrightarrow
x+y=0.
```

于是

```math
\boxed{
N(J)=
\left\{\begin{bmatrix}t\\-t\end{bmatrix}:t\in\mathbb R\right\}
=\mathrm{span}\left\{\begin{bmatrix}1\\-1\end{bmatrix}\right\}
}.
```

这就是二维平面中的直线 $`y=-x`$。例如

```math
J\begin{bmatrix}3\\-3\end{bmatrix}
=\begin{bmatrix}0\\0\end{bmatrix}.
```

$`(-1,1)^T`$ 也可以作为这条线的基，与 $`(1,-1)^T`$ 只差非零倍数。

**列空间记录所有可能的输出，零空间记录哪些输入被送到零。** 本例的两条直线恰好正交，但一般矩阵的列空间与零空间不一定正交，甚至可能处于不同维数的环境空间。

### 3. 特征值是 2 和 0，不是 1 和 0

```math
\det(J-\lambda I_2)
=(1-\lambda)^2-1
=\lambda(\lambda-2).
```

所以

```math
\lambda=2:\quad
E_2=\mathrm{span}\left\{\begin{bmatrix}1\\1\end{bmatrix}\right\},
```

```math
\lambda=0:\quad
E_0=N(J)=\mathrm{span}\left\{\begin{bmatrix}1\\-1\end{bmatrix}\right\}.
```

因为

```math
J\begin{bmatrix}1\\1\end{bmatrix}
=2\begin{bmatrix}1\\1\end{bmatrix},
\qquad
J\begin{bmatrix}1\\-1\end{bmatrix}
=\mathbf0.
```

不是“输出落在一条线上”就一定是投影矩阵；还要看已经位于目标空间的向量是否保持不变。

### 4. 精确的动作分解：先正交投影，再放大 2 倍

令 $`a=(1,1)^T`$。上一章的投影公式给出

```math
P=\frac{aa^T}{a^Ta}
=\frac12\begin{bmatrix}1&1\\1&1\end{bmatrix}.
```

这是投影到 $`y=x`$ 的正交投影矩阵，且

```math
\boxed{J=2P}.
```

对任意 $`u=(x,y)^T`$，

```math
Pu=
\begin{bmatrix}
(x+y)/2\\
(x+y)/2
\end{bmatrix},
\qquad
Ju=2Pu.
```

例如 $`u=(2,0)^T`$：

```math
Pu=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
Ju=\begin{bmatrix}2\\2\end{bmatrix}.
```

因此 $`J`$ 的动作是：先沿垂直于 $`y=x`$ 的方向投到这条线上，再把得到的向量相对原点放大 2 倍。

```math
P^2=P,\qquad J^2=2J\neq J.
```

所以 $`P`$ 是投影矩阵，而 $`J`$ 不是。二者有相同的列空间和零空间，但在保留方向上的缩放不同：

| 矩阵 | 在 $`y=x`$ 上 | 在 $`y=-x`$ 上 | 特征值 |
| --- | --- | --- | --- |
| $`P`$ | 保持不变 | 归零 | 1、0 |
| $`J=2P`$ | 放大 2 倍 | 归零 | 2、0 |

这也呼应前面的结论：**基本子空间相同，不代表矩阵对每个向量的具体作用相同。**

## 六、投影矩阵：怎样从矩阵找出目标空间与垂直方向

### 1. 先区分 projection 与 orthogonal projection

满足 $`P^2=P`$ 的方阵称为投影矩阵，即 **idempotent（幂等）**矩阵。若还是实对称矩阵，

```math
\boxed{P^2=P,\qquad P^T=P},
```

则是 **orthogonal projection matrix（正交投影矩阵）**。

本章中“垂直投到平面”的几何图像专指正交投影。一般的投影也满足下面的 0、1 特征值结论，但其零空间不一定垂直于目标空间。

### 2. 为什么特征值只能取 0 或 1

若 $`Pv=\lambda v`$ 且 $`v\neq\mathbf0`$，由 $`P^2=P`$，

```math
P^2v=P(\lambda v)=\lambda Pv=\lambda^2v,
```

另一方面 $`P^2v=Pv=\lambda v`$。所以

```math
(\lambda^2-\lambda)v=\mathbf0.
```

因为 $`v\neq\mathbf0`$，至少一个分量非零，于是

```math
\boxed{\lambda(\lambda-1)=0
\Longrightarrow\lambda\in\{0,1\}}.
```

“只能取 0 或 1”不表示二者必定同时出现：$`P=I_n`$ 只有特征值 1，$`P=0`$ 只有特征值 0。投影到非零真子空间时才两者都有。

### 3. 为什么目标空间正好等于列空间，也等于 $`N(P-I)`$

令投影目标空间为 $`S=\mathrm{Col}(P)`$。严格证明

```math
\boxed{N(P-I_n)=\mathrm{Col}(P)}.
```

首先，若 $`y\in\mathrm{Col}(P)`$，则存在 $`x`$ 使 $`y=Px`$。于是

```math
Py=P(Px)=P^2x=Px=y.
```

所以 $`y\in N(P-I_n)`$。

反过来，若 $`(P-I_n)y=\mathbf0`$，则 $`y=Py`$，它本身就是一个矩阵 $`P`$ 的输出，所以 $`y\in\mathrm{Col}(P)`$。

因此在 $`1`$ 确为特征值时，

```math
E_1=N(P-I_n)=S.
```

这不只是说“输出落在目标空间”：目标空间中的每个 $`y`$ 都可以用它自己作为输入产生，即 $`Py=y`$，因而所有输出确实覆盖整个目标空间。

### 4. 为什么被消去的空间在正交投影中是垂直方向

零特征值对应的解集是 $`E_0=N(P)`$。对正交投影，

```math
\boxed{N(P)=S^\perp}.
```

这是四个基本子空间的正交关系在 $`P^T=P`$ 时的应用：

```math
N(P)=\mathrm{Col}(P^T)^\perp
=\mathrm{Col}(P)^\perp.
```

任何输入 $`x`$ 都能分解为

```math
\boxed{x=Px+(I_n-P)x}.
```

其中 $`Px\in S`$，且

```math
P(I_n-P)x=(P-P^2)x=\mathbf0,
```

所以 $`(I_n-P)x\in N(P)`$。当 $`P`$ 是正交投影时，这两部分还互相正交。

若一个向量的这两部分都非零，它既不在固定空间，也不在零空间，因此不是投影矩阵的特征向量。不能把所有向量都归入“在目标空间内”或“垂直于目标空间”这两类；一般向量是这两类分量的和。

### 5. 已知矩阵，实际怎样求基

要找目标空间，可以任选一种方法：

1. 解 $`(P-I_n)v=\mathbf0`$，写出零空间的一组基。
2. 求 $`\mathrm{Col}(P)`$ 的一组基：通过消元找出主元列编号，再取**原矩阵 $`P`$ 的对应列**。

必须选出 $`\mathrm{rank}(P)`$ 个线性无关列，足以张成整个目标空间，而不是随意找到一小组独立列就停止。

要找被消去的空间，解 $`Pv=\mathbf0`$，用自由变量写出 $`N(P)`$ 的基。若 $`P`$ 是正交投影，这同时就是目标空间正交补的基。

目标空间维数是 $`r=\mathrm{rank}(P)`$，被消去空间维数是 $`n-r`$。不一定总是“一个平面和一条直线”，那是 $`n=3,r=2`$ 时的特例。

### 6. 三维算例：不与坐标平面重合的投影平面

考虑

```math
P=
\begin{bmatrix}
1/2&1/2&0\\
1/2&1/2&0\\
0&0&1
\end{bmatrix}.
```

可直接验证 $`P^T=P`$、$`P^2=P`$。它作用在 $`(x,y,z)^T`$ 上得到

```math
P\begin{bmatrix}x\\y\\z\end{bmatrix}
=
\begin{bmatrix}(x+y)/2\\(x+y)/2\\z\end{bmatrix}.
```

**求目标平面。** 解 $`(P-I_3)v=\mathbf0`$，只得到约束 $`x=y`$，而 $`z`$ 自由：

```math
S=E_1=
\left\{\begin{bmatrix}s\\s\\t\end{bmatrix}:s,t\in\mathbb R\right\}
=\mathrm{span}\left\{
\begin{bmatrix}1\\1\\0\end{bmatrix},
\begin{bmatrix}0\\0\\1\end{bmatrix}
\right\}.
```

也可以直接看原矩阵：第二列与第一列相同，第一列和第三列线性无关，二者张成同一个平面。第一列乘 2 后就是上式选取的第一个基向量。

**求垂直方向。** 解 $`Pv=\mathbf0`$，得到 $`x+y=0,z=0`$：

```math
E_0=N(P)=
\mathrm{span}\left\{\begin{bmatrix}1\\-1\\0\end{bmatrix}\right\}.
```

检查正交：

```math
\begin{bmatrix}1\\-1\\0\end{bmatrix}^{T}
\begin{bmatrix}s\\s\\t\end{bmatrix}=s-s=0.
```

所以三维输入空间被分成固定平面 $`x=y`$ 与它的法线：

```math
\mathbb R^3=E_1\oplus E_0,\qquad
\dim E_1+\dim E_0=2+1=3.
```

若换成更简单的 $`P=\mathrm{diag}(1,1,0)`$，同样过程给出目标平面 $`z=0`$、垂直方向 $`\mathrm{span}\{(0,0,1)^T\}`$。

## 七、特征值之和与积：怎样连接 trace 和 determinant

### 1. 先约定特征多项式的符号与重数

求解时常用 $`\det(A-\lambda I_n)=0`$。为了比较多项式系数，这里采用首项系数为 1 的约定：

```math
\boxed{p_A(\lambda)=\det(\lambda I_n-A)}.
```

因为

```math
\det(A-\lambda I_n)=(-1)^n p_A(\lambda),
```

两种写法的根完全相同。

若允许复数根，$`n`$ 次特征多项式可以分解为

```math
p_A(\lambda)=\prod_{j=1}^n(\lambda-\lambda_j).
```

这里列出 $`n`$ 个根时，**按 algebraic multiplicity（代数重数）重复计数**，不是说一定有 $`n`$ 个不同特征值或 $`n`$ 个线性无关特征向量。

例如 $`A=2I_2`$ 时，$`p_A(\lambda)=(\lambda-2)^2`$，求和、求积都要把 2 计入两次。特征空间的维数则叫 geometric multiplicity（几何重数），是另一个概念。

### 2. 二阶：用熟悉的韦达定理完整推导

设

```math
A=\begin{bmatrix}a&b\\c&d\end{bmatrix},
\qquad
\mathrm{tr}(A)=a+d.
```

trace（迹）就是主对角线元素的和。展开：

```math
\begin{aligned}
p_A(\lambda)
&=(\lambda-a)(\lambda-d)-bc\\
&=\lambda^2-(a+d)\lambda+(ad-bc)\\
&=\lambda^2-\mathrm{tr}(A)\lambda+\det(A).
\end{aligned}
```

又因为它的根是 $`\lambda_1,\lambda_2`$，

```math
p_A(\lambda)
=(\lambda-\lambda_1)(\lambda-\lambda_2)
=\lambda^2-(\lambda_1+\lambda_2)\lambda+\lambda_1\lambda_2.
```

逐项比较，得到

```math
\boxed{
\lambda_1+\lambda_2=\mathrm{tr}(A),\qquad
\lambda_1\lambda_2=\det(A)
}.
```

这是“根与系数”的结论，不需要假设矩阵已存在一组特征向量基。

### 3. 高阶：次高次项为何只看主对角线

对 $`n\geq2`$，展开 $`\det(\lambda I_n-A)`$。

主对角线乘积为

```math
\prod_{i=1}^n(\lambda-a_{ii}).
```

取每个因子中的 $`\lambda`$，产生 $`\lambda^n`$；取一个 $`-a_{ii}`$、其余取 $`\lambda`$，产生

```math
-\left(\sum_{i=1}^n a_{ii}\right)\lambda^{n-1}.
```

为什么其他 determinant 展开项不会贡献 $`\lambda^{n-1}`$？因为每项要从每行、每列各取一个元素；一旦不沿主对角线选取，就至少有两处离开主对角线。非对角元素不含 $`\lambda`$，所以这些项的次数最多为 $`n-2`$。

因此 $`\lambda^{n-1}`$ 的系数恰好是 $`-\mathrm{tr}(A)`$。另一方面，由根的乘积分解，该系数为 $`-\sum_j\lambda_j`$。比较得

```math
\boxed{\sum_{j=1}^n\lambda_j=\mathrm{tr}(A)}.
```

再比较常数项：

```math
p_A(0)=\det(-A)=(-1)^n\det(A),
```

```math
\prod_{j=1}^n(0-\lambda_j)
=(-1)^n\prod_{j=1}^n\lambda_j.
```

所以

```math
\boxed{\prod_{j=1}^n\lambda_j=\det(A)}.
```

一阶情形直接由 $`p_A(\lambda)=\lambda-a_{11}`$ 得到同样结论。

注意符号：若使用 $`\det(A-\lambda I_n)`$，首项是 $`(-1)^n\lambda^n`$，次高次项系数是 $`(-1)^{n-1}\mathrm{tr}(A)`$，不能不看奇偶性就说“该系数等于 trace”。

### 4. 用于检查答案，但不能替代完整求解

第四节矩阵特征值为 3、1：

```math
3+1=4=\mathrm{tr}(A),\qquad
3\cdot1=3=\det(A).
```

全 1 矩阵 $`J`$ 特征值为 2、0：

```math
2+0=2=\mathrm{tr}(J),\qquad
2\cdot0=0=\det(J).
```

一个特征值描述一个特征方向上的缩放，determinant 描述有向体积缩放；它们一般不是同一个数，但也不是绝对不可能数值相等。

特别地：

```math
\boxed{
0\text{ 是 }A\text{ 的特征值}
\Longleftrightarrow
N(A)\neq\{\mathbf0\}
\Longleftrightarrow
A\text{ singular}
}.
```

求和与求积是有力的检查工具，但尤其在高维时，不能仅凭这两个数确定所有特征值，更不能确定全部特征向量。

### 5. 为什么写出“两个 2”：重根不等于两个独立方向

课堂例子是

```math
T=\begin{bmatrix}2&1\\0&2\end{bmatrix}.
```

用二阶行列式公式：

```math
\begin{aligned}
\det(T-\lambda I_2)
&=\det\begin{bmatrix}2-\lambda&1\\0&2-\lambda\end{bmatrix}\\
&=(2-\lambda)(2-\lambda)-1\cdot0\\
&=(2-\lambda)^2.
\end{aligned}
```

所以特征方程只有一个不同的解 $`\lambda=2`$，但这个根的代数重数为 2，按重数排列就写成

```math
\boxed{\lambda_1=2,\qquad\lambda_2=2}.
```

“两个 2”是在记录特征多项式的两个根位置，不是说出现两个不同数值，也没有自动提供两个独立特征向量。要知道有多少独立方向，还必须解

```math
(T-2I_2)v=
\begin{bmatrix}0&1\\0&0\end{bmatrix}
\begin{bmatrix}a\\b\end{bmatrix}
=\mathbf0.
```

这只要求 $`b=0`$，$`a`$ 自由，因此

```math
E_2(T)=
\mathrm{span}\left\{\begin{bmatrix}1\\0\end{bmatrix}\right\},
\qquad
\dim E_2(T)=1.
```

所有非零解彼此成倍数，无法选出两根线性无关的特征向量。把同一个向量重复放两列组成的 $`S`$ 不可逆，所以 $`T`$ 不能对角化。

与前面出现的 $`2I_2`$ 对照：

| 矩阵 | 特征多项式 | 不同特征值 | 对应特征空间维数 | 能否对角化 |
| --- | --- | --- | --- | --- |
| $`2I_2`$ | $`(\lambda-2)^2`$ | 只有 2 | 2 | 能，原本就是对角矩阵 |
| $`T`$ | $`(\lambda-2)^2`$ | 只有 2 | 1 | 不能 |

因此**有重根不一定不能对角化；本例不能对角化，是因为独立特征向量不足。** 同时 $`\mathrm{tr}(T)=2+2=4`$、$`\det(T)=2\cdot2=4`$，再次说明求和与求积必须计入重数。

### 6. 实矩阵不一定有实特征方向

平面旋转 90 度的矩阵

```math
R=\begin{bmatrix}0&-1\\1&0\end{bmatrix}
```

满足

```math
p_R(\lambda)=\lambda^2+1.
```

在实数域没有根，也没有非零实特征向量；在复数域有特征值 $`i,-i`$。仍有

```math
i+(-i)=0=\mathrm{tr}(R),\qquad
i(-i)=1=\det(R).
```

因此“$`n\times n`$ 矩阵有 $`n`$ 个特征值”必须说明：是在复数域中，按代数重数计数，而不是一定有 $`n`$ 条实特征直线。

## 八、为什么研究特征值：重复变换与坐标简化

### 1. 在一个特征方向上，矩阵幂变成标量幂

若 $`Av=\lambda v`$，则

```math
A^2v=A(\lambda v)=\lambda Av=\lambda^2v.
```

依次得到，对正整数 $`k`$，

```math
\boxed{A^kv=\lambda^kv}.
```

对实特征方向而言：

| 条件 | 该方向上反复变换的效果 |
| --- | --- |
| $`\lvert\lambda\rvert<1`$ | 衰减到零 |
| $`\lvert\lambda\rvert>1`$ | 长度增长 |
| $`\lambda=1`$ | 保持不变 |
| $`\lambda=-1`$ | 在正反方向间交替 |

这些首先是对该特征方向的结论。不能只看到一个 $`\lambda=1`$，就断定所有初始向量都保持不变。

### 2. 全 1 矩阵怎样分开处理两个分量

任意二维向量都能写成

```math
\begin{bmatrix}x\\y\end{bmatrix}
=
\frac{x+y}{2}\begin{bmatrix}1\\1\end{bmatrix}
+
\frac{x-y}{2}\begin{bmatrix}1\\-1\end{bmatrix}.
```

这两部分分别落在 $`J`$ 的 $`E_2`$ 与 $`E_0`$ 中。于是 $`k\geq1`$ 时，

```math
J^k\begin{bmatrix}x\\y\end{bmatrix}
=
2^k\frac{x+y}{2}\begin{bmatrix}1\\1\end{bmatrix}.
```

另一分量第一次就被消去，所以也可以写成

```math
\boxed{J^k=2^{k-1}J\qquad(k\geq1)}.
```

若输入本来就在 $`N(J)`$ 中，它会直接归零，而不是因为矩阵存在特征值 2 就一定增长。**是否激发某个特征方向，还取决于初始向量在该方向上的分量。**

### 3. 特征向量矩阵 $`S`$：为什么初始向量可以写成 $`Sc`$

假设 $`A`$ 有 $`n`$ 个线性无关的特征向量 $`v_1,\ldots,v_n`$，它们构成所讨论空间的一组基。板书把这些向量记为 $`x_i`$；本文统一用 $`v_i`$，避免和向量的分量混淆。

把它们按列排列，并把相应特征值按**同一顺序**放在对角线上：

```math
S=[v_1\ \cdots\ v_n],\qquad
\Lambda=\mathrm{diag}(\lambda_1,\ldots,\lambda_n).
```

任意初始向量都能唯一写成

```math
u_0=c_1v_1+\cdots+c_nv_n.
```

令

```math
c=\begin{bmatrix}c_1\\\vdots\\c_n\end{bmatrix},
```

按照矩阵乘向量的列组合含义，

```math
\boxed{u_0=Sc,\qquad c=S^{-1}u_0}.
```

这里没有新的运算，只是把“方向的线性组合”写成“矩阵乘系数向量”。

| 符号 | 含义 | 尺寸 |
| --- | --- | --- |
| $`u_0`$ | 已知的初始状态，采用原坐标 | $`n\times1`$ |
| $`v_i`$ | 特征向量基中的一个向量 | $`n\times1`$ |
| $`S`$ | 特征向量按列组成的换基矩阵 | $`n\times n`$ |
| $`c`$ | 初始状态在特征向量基下的坐标 | $`n\times1`$ |
| $`\Lambda`$ | 与 $`S`$ 各列依次配对的特征值 | $`n\times n`$ |

因此 $`u_0`$ 和 $`c`$ 不是两个不同状态，而是同一个状态在不同基下的坐标；一般也不能把 $`c_i`$ 直接当成 $`u_0`$ 原来的第 $`i`$ 个分量。

如果特征向量没有构成整个空间的基，某些特殊的 $`u_0`$ 仍可能写成它们的组合，但不能保证**任意**初始向量都能写出，也不能直接使用 $`S^{-1}`$。

### 4. 对角化与矩阵幂：为什么中间的逆矩阵会消掉

逐列使用 $`Av_i=\lambda_iv_i`$：

```math
AS=[Av_1\ \cdots\ Av_n]
=[\lambda_1v_1\ \cdots\ \lambda_nv_n]
=S\Lambda.
```

由于 $`S`$ 可逆，在右边乘 $`S^{-1}`$：

```math
\boxed{A=S\Lambda S^{-1}}.
```

这就是 **diagonalization（对角化）**。再计算平方：

```math
\begin{aligned}
A^2
&=(S\Lambda S^{-1})(S\Lambda S^{-1})\\
&=S\Lambda(S^{-1}S)\Lambda S^{-1}\\
&=S\Lambda^2S^{-1}.
\end{aligned}
```

同理，对非负整数 $`k`$ 有

```math
\boxed{A^k=S\Lambda^kS^{-1}}.
```

$`k=0`$ 时两边都是单位矩阵；正整数 $`k`$ 时，

```math
\Lambda^k=\mathrm{diag}(\lambda_1^k,\ldots,\lambda_n^k).
```

所以计算路线是：$`S^{-1}`$ 把原坐标换成特征坐标，$`\Lambda^k`$ 分别缩放各个分量，$`S`$ 再组合回原坐标。没有使用矩阵乘法的交换律，也不是对 $`A`$ 的每个元素分别取 $`k`$ 次方。

第七节的 $`T`$ 没有足够的独立特征向量，不能使用这一对角化。另一个前面讨论过的例子

```math
H=\begin{bmatrix}1&1\\0&1\end{bmatrix}
```

也只有特征方向 $`(1,0)^T`$，但

```math
H^k=\begin{bmatrix}1&k\\0&1\end{bmatrix}.
```

因此即使特征值全为 1，某些输入仍随迭代增长；不能把“特征方向上的行为”直接推广成任意输入的行为。

### 5. 从每步递推到一次写出全部状态

若每一步都按同一个矩阵更新，

```math
\boxed{u_{k+1}=Au_k,\qquad k=0,1,2,\ldots},
```

则 $`u_k`$ 是第 $`k`$ 步的整个状态向量，不是单个特征向量，也不是向量的第 $`k`$ 个分量。

连续代入：

```math
u_1=Au_0,\qquad
u_2=Au_1=A^2u_0,\qquad
u_k=A^ku_0.
```

若已经求得 $`u_0=\sum_i c_iv_i`$，利用矩阵乘法的线性性，对 $`k\geq1`$，

```math
\begin{aligned}
u_k
&=A^k\left(\sum_{i=1}^n c_iv_i\right)\\
&=\sum_{i=1}^n c_iA^kv_i\\
&=\boxed{\sum_{i=1}^n c_i\lambda_i^kv_i}.
\end{aligned}
```

原来的 $`c_i`$ 是初始状态的固定系数；随时间变化的特征坐标是 $`c_i\lambda_i^k`$。不是每一步都重新求一遍系数，也不是说每一步的特征坐标一直等于初始 $`c_i`$。

同一个结果也可以写成

```math
\boxed{u_k=S\Lambda^kc=S\Lambda^kS^{-1}u_0}.
```

> **先找特征方向，再按初始状态算权重；时间推进只让每个方向的权重多乘一次自己的特征值。**

### 6. 完整例子：如何计算 $`A^{100}u_0`$

继续使用第四节的矩阵：

```math
A=\begin{bmatrix}2&1\\1&2\end{bmatrix},
\qquad
v_1=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
v_2=\begin{bmatrix}1\\-1\end{bmatrix}.
```

对应 $`\lambda_1=3,\lambda_2=1`$。现在另给初始状态

```math
u_0=\begin{bmatrix}3\\1\end{bmatrix}.
```

**先求初始权重。** 令 $`u_0=c_1v_1+c_2v_2`$，按分量相等列式：

```math
\begin{cases}
c_1+c_2=3,\\
c_1-c_2=1.
\end{cases}
\qquad\Longrightarrow\qquad
c_1=2,\quad c_2=1.
```

也就是

```math
u_0=2v_1+v_2.
```

用矩阵坐标写，完全等价于

```math
S=\begin{bmatrix}1&1\\1&-1\end{bmatrix},
\qquad
S^{-1}=\frac12\begin{bmatrix}1&1\\1&-1\end{bmatrix},
\qquad
c=S^{-1}u_0=\begin{bmatrix}2\\1\end{bmatrix}.
```

**再推进到第 $`k`$ 步。**

```math
u_k=2\cdot3^kv_1+1^kv_2
=\begin{bmatrix}2\cdot3^k+1\\2\cdot3^k-1\end{bmatrix}.
```

因此

```math
\boxed{
A^{100}u_0=
\begin{bmatrix}2\cdot3^{100}+1\\2\cdot3^{100}-1\end{bmatrix}
}.
```

这里求的是 $`A^{100}u_0`$ 这个向量，不必先写出整个 $`A^{100}`$。如需完整矩阵幂，对角化还给出

```math
A^k
=\frac12
\begin{bmatrix}
3^k+1&3^k-1\\
3^k-1&3^k+1
\end{bmatrix}.
```

检查 $`k=0`$ 得 $`I_2`$；$`k=1`$ 时 $`u_1=(7,5)^T`$，直接计算 $`Au_0`$ 也得到同一结果。

这个例子完成了“先拆初始向量，再对每个特征方向求幂”的整条计算链。下一节把同一方法用于一个标量数列：先把相邻两项装进状态向量，就能处理二阶递推。

## 九、斐波那契：从二阶递推到显式通式

这一例沿用上一节的方法。关键问题依次是：怎样选 $`u_k`$，怎样找特征对，怎样拆初值，最后从哪个分量读出答案。

### 1. 为什么把相邻两项装进一个状态向量

约定

```math
F_0=0,\qquad F_1=1,\qquad
F_{k+2}=F_{k+1}+F_k\quad(k\geq0).
```

单独一个 $`F_k`$ 一般不足以决定下一项，还需要相邻的另一项。于是记录完整的更新信息：

```math
\boxed{u_k=\begin{bmatrix}F_{k+1}\\F_k\end{bmatrix}}.
```

那么

```math
\begin{aligned}
u_{k+1}
&=\begin{bmatrix}F_{k+2}\\F_{k+1}\end{bmatrix}
=\begin{bmatrix}F_{k+1}+F_k\\F_{k+1}\end{bmatrix}\\
&=\underbrace{\begin{bmatrix}1&1\\1&0\end{bmatrix}}_{B}
\begin{bmatrix}F_{k+1}\\F_k\end{bmatrix}.
\end{aligned}
```

所以

```math
\boxed{
u_{k+1}=Bu_k,\qquad
u_0=\begin{bmatrix}F_1\\F_0\end{bmatrix}
=\begin{bmatrix}1\\0\end{bmatrix},\qquad
u_k=B^ku_0
}.
```

板书记这个矩阵为 $`A`$；为与上一节的矩阵区分，本文记为 $`B`$。它的第一行“相加”，第二行“把原来的第一项移到第二项位置”，恰好编码数列的更新规则。

### 2. 为什么再次求 $`\det(B-\lambda I)=0`$

要把 $`B^k`$ 变成特征值的乘方，先为这个新矩阵求特征对：

```math
B-\lambda I_2=
\begin{bmatrix}1-\lambda&1\\1&-\lambda\end{bmatrix}.
```

```math
\det(B-\lambda I_2)
=(1-\lambda)(-\lambda)-1
=\lambda^2-\lambda-1.
```

要求这个 determinant 为零，绝不是让 $`B-\lambda I_2`$ 成为零矩阵；它的两个非对角元素一直都是 1。

解二次方程：

```math
\boxed{
\lambda_1=\phi=\frac{1+\sqrt5}{2},\qquad
\lambda_2=\psi=\frac{1-\sqrt5}{2}
}.
```

后面会用到：

```math
\phi+\psi=1,\qquad
\phi\psi=-1,\qquad
\phi-\psi=\sqrt5,
```

```math
\phi^2=\phi+1,\qquad\psi^2=\psi+1.
```

前两个也与 $`\mathrm{tr}(B)=1`$、$`\det(B)=-1`$ 一致。

### 3. 已知特征值后，怎样算出特征向量

对每个已经求出的 $`\lambda`$，写 $`v=(a,b)^T`$ 并解

```math
\begin{bmatrix}1-\lambda&1\\1&-\lambda\end{bmatrix}
\begin{bmatrix}a\\b\end{bmatrix}=\mathbf0.
```

按分量展开：

```math
\begin{cases}
(1-\lambda)a+b=0,\\
a-\lambda b=0.
\end{cases}
```

第二式给出 $`a=\lambda b`$。如果 $`b=0`$，则 $`a=0`$，不能成为特征向量；因此可以用非零缩放把 $`b`$ 取成 1，选

```math
v=\begin{bmatrix}\lambda\\1\end{bmatrix}.
```

它还满足第一式，因为

```math
(1-\lambda)\lambda+1
=-(\lambda^2-\lambda-1)=0.
```

这里正是使用了 $`\lambda`$ 已满足特征方程。一般解为 $`t(\lambda,1)^T`$，取 $`t\neq0`$ 得到对应特征向量。因此可选

```math
\boxed{
v_1=\begin{bmatrix}\phi\\1\end{bmatrix},\qquad
v_2=\begin{bmatrix}\psi\\1\end{bmatrix}
}.
```

也可直接核对：

```math
B\begin{bmatrix}\lambda\\1\end{bmatrix}
=\begin{bmatrix}\lambda+1\\\lambda\end{bmatrix}
=\begin{bmatrix}\lambda^2\\\lambda\end{bmatrix}
=\lambda\begin{bmatrix}\lambda\\1\end{bmatrix}.
```

把它们放入

```math
S=\begin{bmatrix}\phi&\psi\\1&1\end{bmatrix},
```

有 $`\det(S)=\phi-\psi=\sqrt5\neq0`$，所以两根向量确实线性无关，足以构成二维空间的一组基。

### 4. 特征向量怎样与 $`u_0`$ 结合：求初始权重

找到特征对后，还需要把指定的初值拆开：

```math
\begin{bmatrix}1\\0\end{bmatrix}
=c_1\begin{bmatrix}\phi\\1\end{bmatrix}
+c_2\begin{bmatrix}\psi\\1\end{bmatrix}.
```

上下分量分别相等：

```math
\begin{cases}
c_1\phi+c_2\psi=1,\\
c_1+c_2=0.
\end{cases}
```

由 $`c_2=-c_1`$，得到 $`c_1(\phi-\psi)=1`$。因此

```math
\boxed{c_1=\frac1{\sqrt5},\qquad c_2=-\frac1{\sqrt5}}.
```

也就是

```math
u_0=\frac1{\sqrt5}v_1-\frac1{\sqrt5}v_2.
```

用坐标变换写：

```math
S^{-1}=\frac1{\sqrt5}
\begin{bmatrix}1&-\psi\\-1&\phi\end{bmatrix},
\qquad
c=S^{-1}u_0=\frac1{\sqrt5}\begin{bmatrix}1\\-1\end{bmatrix}.
```

系数由初值决定；如果数列的前两项改变，递推矩阵和特征对可以不变，但 $`c_1,c_2`$ 必须重新计算。

### 5. 先求整个状态，再读出斐波那契通式

利用 $`B^kv_i=\lambda_i^kv_i`$：

```math
\begin{aligned}
u_k
&=B^ku_0\\
&=\frac{\phi^k}{\sqrt5}v_1-\frac{\psi^k}{\sqrt5}v_2\\
&=\frac1{\sqrt5}
\left(
\phi^k\begin{bmatrix}\phi\\1\end{bmatrix}
-\psi^k\begin{bmatrix}\psi\\1\end{bmatrix}
\right)\\
&=\begin{bmatrix}
\dfrac{\phi^{k+1}-\psi^{k+1}}{\sqrt5}\\[6pt]
\dfrac{\phi^k-\psi^k}{\sqrt5}
\end{bmatrix}.
\end{aligned}
```

一开始定义的是 $`u_k=(F_{k+1},F_k)^T`$，所以**第二个分量才是 $`F_k`$，第一个分量是 $`F_{k+1}`$**。逐项比较得到 **Binet's formula（比内公式）**：

```math
\boxed{
F_k=\frac{\phi^k-\psi^k}{\sqrt5}
=\frac1{\sqrt5}
\left[
\left(\frac{1+\sqrt5}{2}\right)^k
-\left(\frac{1-\sqrt5}{2}\right)^k
\right],\qquad k\geq0
}.
```

若改用 $`(F_k,F_{k-1})^T`$ 作为状态，也能建立递推，但初值和步数编号会相应改变。不能在同一推导中混用两种下标约定。

### 6. 检查初值、递推与整个矩阵幂

**检查初值和下一项。**

```math
F_0=\frac{1-1}{\sqrt5}=0,\qquad
F_1=\frac{\phi-\psi}{\sqrt5}=1,
```

```math
F_2=\frac{(\phi-\psi)(\phi+\psi)}{\sqrt5}=1.
```

**检查递推。** 因为 $`\phi^2=\phi+1`$、$`\psi^2=\psi+1`$，分别乘相应的 $`k`$ 次方再相减：

```math
\frac{\phi^{k+2}-\psi^{k+2}}{\sqrt5}
=
\frac{\phi^{k+1}-\psi^{k+1}}{\sqrt5}
+
\frac{\phi^k-\psi^k}{\sqrt5}.
```

这就是原数列的递推。相同初值与相同递推保证两者逐项一致；表达式含根号不妨碍结果为整数。

| $`k`$ | $`u_k=(F_{k+1},F_k)^T`$ | 读出的 $`F_k`$ |
| --- | --- | --- |
| 0 | $`(1,0)^T`$ | 0 |
| 1 | $`(1,1)^T`$ | 1 |
| 2 | $`(2,1)^T`$ | 1 |
| 3 | $`(3,2)^T`$ | 2 |
| 4 | $`(5,3)^T`$ | 3 |
| 5 | $`(8,5)^T`$ | 5 |

还可以从状态向量得到整个矩阵幂。由于 $`u_0=e_1`$，$`u_k=B^ke_1`$ 正是 $`B^k`$ 的第一列；又因为 $`Be_2=e_1`$，当 $`k\geq1`$，

```math
B^ke_2=B^{k-1}e_1=u_{k-1}.
```

所以

```math
\boxed{
B^k=\begin{bmatrix}
F_{k+1}&F_k\\
F_k&F_{k-1}
\end{bmatrix}\qquad(k\geq1)
}.
```

$`k=0`$ 时使用 $`B^0=I_2`$，不要在未定义负下标数列时直接代入 $`F_{-1}`$。

### 7. 这一例最终解释了什么

两个特征值满足 $`\phi>1`$、$`|\psi|<1`$。本题初值在 $`\phi`$ 对应方向上的系数非零，所以长期增长由 $`\phi^k`$ 主导，$`\psi^k`$ 则衰减。

整条任务链是：

```text
F_(k+2)=F_(k+1)+F_k
    -> 两项决定下一项，选 u_k=(F_(k+1),F_k)^T
    -> u_(k+1)=B u_k，初值 u_0=(1,0)^T
    -> 解 det(B-λI)=0，得到 φ、ψ
    -> 分别解 (B-λI)v=0，得到 v_1、v_2
    -> 解 S c=u_0，得到初始权重
    -> u_k=S Λ^k c，把矩阵迭代变成数的乘方
    -> 读取第二个分量，得到 F_k
```

> **特征值决定各方向怎样随时间变化；特征向量决定这些方向是什么；初值决定各方向占多少；状态定义决定最后从哪里读出所求数列。**

## 十、常见混淆与修正

| 容易混淆的说法 | 应怎样理解 |
| --- | --- |
| “$`Av=\lambda v`$ 已知 $`\lambda`$，求 $`v`$” | 一般已知 $`A`$，同时找 $`\lambda,v`$；若已给 $`\lambda`$，才直接求移位零空间 |
| “为了求解，人为要求行列式等于零” | 这是存在非零特征向量的充要条件 |
| “有特征值说明 $`A`$ 不可逆” | 必须是 $`0`$ 为特征值才说明 $`A`$ 不可逆；一般是 $`A-\lambda I`$ 不可逆 |
| “singular 与 defective 是同一个词” | singular 是不可逆；defective 指缺少足够独立特征向量，不是同一概念 |
| “特征向量的任何倍数仍是特征向量” | 必须是非零倍数；零向量被定义排除 |
| “每个特征空间只有一条方向” | 特征空间可以是直线、平面，甚至整个空间 |
| “矩阵的四个元素表示四维坐标” | $`2\times2`$ 矩阵在本章表示二维到二维的变换，其两列是两个标准基的输出 |
| “一个向量就是一条直线” | 一个向量给出一根箭头，它的所有实数倍构成的 span 才是直线 |
| “向量 $`v`$ 在两边相同，所以可以约掉” | 单个输入上作用相同，不等于两个矩阵在所有输入上相同 |
| “所有投影都沿垂直方向” | 正交投影还需对称性 $`P^T=P`$；一般投影只要求 $`P^2=P`$ |
| “投影矩阵必定同时有 0、1 两个特征值” | 投到非零真子空间时如此；零矩阵和单位矩阵是边界情况 |
| “像全 1 矩阵那样输出到直线，就是投影” | $`J=2P`$ 还沿保留方向放大，$`J^2\neq J`$，不是投影矩阵 |
| “$`\lambda`$ 就是 $`\det(A)`$” | determinant 是全部特征值按代数重数的乘积，不是一般的单个特征值 |
| “特征值之和等于 trace，不用考虑重数和复数” | 要统计特征多项式的全部复根并计入代数重数 |
| “找到特征值就一定能用特征向量作一组基” | 还要检查是否有 $`n`$ 个线性无关特征向量 |
| “$`\lambda_1=\lambda_2=2`$ 表示两个不同特征值或两根独立特征向量” | 一个不同特征值按代数重数计两次；独立方向要由 $`N(A-2I)`$ 求出 |
| “有重复特征值就一定不能对角化” | $`2I_2`$ 有重根且能对角化；还需检查独立特征向量的总数 |
| “求特征值要使 $`A-\lambda I=0`$” | 是 $`\det(A-\lambda I)=0`$，不是矩阵本身为零 |
| “$`u_0`$ 的分量就是 $`c_i`$” | $`c=S^{-1}u_0`$ 是特征基下的坐标，一般不等于原坐标 |
| “已知特征向量，递推题就做完了” | 还需按初值求系数，再按特征值求幂 |
| “$`A^k`$ 是矩阵每个元素分别取 $`k`$ 次方” | 是矩阵连续相乘 $`k`$ 次；可对角化时改写为 $`S\Lambda^kS^{-1}`$ |
| “$`u_k`$ 的第一个分量总是 $`F_k`$” | 由状态定义决定；本节 $`u_k=(F_{k+1},F_k)^T`$，第二分量才是 $`F_k`$ |

## 十一、自测与本章总结

### 1. 七个检查理解的问题

**题一。** $`A=\mathrm{diag}(2,3)`$，取 $`v=(1,0)^T`$，能从 $`Av=2v`$ 推出 $`A=2I_2`$ 吗？

不能。该等式只约束一个方向；$`A`$ 在第二个标准方向上乘 3。

**题二。** $`A=2I_2`$ 时，$`\lambda=2`$ 的特征空间是一条线吗？

不是，$`A-2I_2=0`$，所以 $`E_2=\mathbb R^2`$。全部非零二维向量都是对应特征向量。

**题三。** $`J=\begin{bmatrix}1&1\\1&1\end{bmatrix}`$ 是投影矩阵吗？它的零空间是哪条线？

不是，因为 $`J^2=2J\neq J`$。零空间为 $`\mathrm{span}\{(1,-1)^T\}`$，即 $`y=-x`$；$`J/2`$ 才是到 $`y=x`$ 的正交投影矩阵。

**题四。** 对第六节的三维投影矩阵，$`u=(2,0,3)^T`$ 是否是特征向量？

```math
Pu=\begin{bmatrix}1\\1\\3\end{bmatrix}.
```

不是。若它等于 $`\lambda u`$，第二个分量必须为 0，与实际输出的 1 矛盾。它的分解是

```math
u=
\begin{bmatrix}1\\1\\3\end{bmatrix}
+
\begin{bmatrix}1\\-1\\0\end{bmatrix},
```

两部分分别属于 $`E_1`$ 和 $`E_0`$。

**题五。** $`H=\begin{bmatrix}1&1\\0&1\end{bmatrix}`$ 的两个特征根都是 1，能否因此推出 $`H=I_2`$？

不能。它只有一个独立特征方向，且 $`H(0,1)^T=(1,1)^T`$。特征值相同不代表矩阵相同，更不代表它在每个向量上都作用为恒等变换。

**题六。** $`T=\begin{bmatrix}2&1\\0&2\end{bmatrix}`$ 与 $`2I_2`$ 都把特征根记为 2、2，为什么前者不能对角化、后者能？

因为 $`\dim N(T-2I_2)=1`$，而 $`\dim N(2I_2-2I_2)=2`$。重数相同，不代表独立特征方向数量相同。

**题七。** 斐波那契状态 $`u_k=(F_{k+1},F_k)^T`$ 满足 $`u_0=(1,0)^T`$。若 $`B^5u_0=(8,5)^T`$，$`F_5`$ 是 8 还是 5？

是 5。第一分量是 $`F_6=8`$，第二分量才是 $`F_5=5`$。读答案必须遵守最初的状态定义。

### 2. 最值得保留的四条主线

**求解主线：**

```math
\boxed{
Av=\lambda v,\ v\neq\mathbf0
\Longleftrightarrow
(A-\lambda I_n)v=\mathbf0,\ v\neq\mathbf0
}.
```

先用 $`\det(A-\lambda I_n)=0`$ 求特征值，再用 $`E_\lambda=N(A-\lambda I_n)`$ 求整个特征空间。

**子空间主线：**

```math
E_0=N(A)\quad\text{当 }0\text{ 是特征值时};
```

对正交投影 $`P`$，

```math
N(P-I_n)=\mathrm{Col}(P),\qquad
N(P)=\mathrm{Col}(P)^\perp.
```

这把新概念重新接回零空间、列空间、基、维数和正交补。

**矩阵结构与计算主线：**

```math
\boxed{
\sum_{j=1}^n\lambda_j=\mathrm{tr}(A),\qquad
\prod_{j=1}^n\lambda_j=\det(A)
}
```

按代数重数计入全部复特征根；在每个特征向量上，还有 $`A^kv=\lambda^kv`$。若进一步有一组特征向量基，才可以据此建立整个空间的特征坐标系并进行对角化。

**递推应用主线：**

```math
u_{k+1}=Au_k
\Longrightarrow
u_k=A^ku_0
=S\Lambda^kS^{-1}u_0.
```

最后一个等号以 $`A`$ 有一组特征向量基为前提。实际操作是“选状态、求特征对、拆初值、求幂、读分量”，不是只求完特征值就结束。

> **本章的终点：寻找的不是“一个神秘数字”，而是矩阵在哪些非零输入方向上可以被一个数替代。特征方程找到这些数，移位零空间给出这些方向；有足够独立方向时，再由初值决定权重，把矩阵迭代转成各特征坐标的独立演化。**

