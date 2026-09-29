# Determinants

## 本章主线与范围

前面已经用 rank、pivot、null space 和列空间判断方阵是否可逆。本章引入一个新的标量：**determinant（行列式）**。它把这些结构压缩成一个数，并同时描述线性变换对有向面积或体积的缩放。

```text
det(A) 的基本性质
    -> 推出相同行、倍数行和零行时 det(A)=0
    -> 理解三种行操作怎样改变 determinant
    -> 用高斯消元把 A 化为三角矩阵 U
    -> determinant = 换行符号 × pivot 乘积
    -> det(A)=0 恰好对应 singular
    -> LU 把消元过程记录成矩阵分解
    -> 推出 det(AB)、det(A^T)、det(cA) 等性质
```

本章与上一章的两种分解要区分：

| 分解 | 来源 | 主要作用 |
| --- | --- | --- |
| $`A=QR`$ | Gram–Schmidt 正交化 | 用标准正交基描述列空间 |
| $`A=LU`$ | Gaussian elimination 高斯消元 | 记录消元倍数并快速求解方程或 determinant |

> **先记住：$`\det(A)`$ 是一个数，不是一个矩阵；通常只对方阵定义 determinant。**

## 一、determinant 是什么

### 1. 输入是方阵，输出是标量

若

```math
A\in\mathbb{R}^{n\times n},
```

则

```math
\det(A)\in\mathbb{R}.
```

因此 $`\det`$ 是一个把方阵映射到实数的函数：

```math
\det:\mathbb{R}^{n\times n}\longrightarrow\mathbb{R}.
```

普通的 $`m\times n`$ 矩形矩阵没有通常意义下的 determinant。它可以讨论 rank、null space 或列空间，但不能直接写 $`\det(A)`$。

### 2. 二阶行列式

对于

```math
A=
\begin{bmatrix}
a&b\\
c&d
\end{bmatrix},
```

有

```math
\boxed{
\det(A)=ad-bc
}.
```

例如

```math
A=
\begin{bmatrix}
2&1\\
4&3
\end{bmatrix},
```

则

```math
\det(A)=2\cdot3-1\cdot4=2.
```

### 3. 几何意义：有向体积缩放因子

若 $`A`$ 是二维线性变换，则 $`|\det(A)|`$ 表示面积被放大的倍数；在三维中，它表示体积被放大的倍数。

- $`|\det(A)|>1`$：面积或体积被放大；
- $`0<|\det(A)|<1`$：面积或体积被压缩；
- $`\det(A)<0`$：除缩放外还发生了方向翻转；
- $`\det(A)=0`$：某个维度被压扁，体积塌为零。

最后一种情况正对应矩阵丢失信息、不可逆。

## 二、determinant 的三个基本性质

在课程中，可以从下面三条基本性质出发理解 determinant。后面很多规则都不是新的公理，而是由这些性质推导出来的。

### 1. 单位矩阵的 determinant 是 1

```math
\boxed{\det(I_n)=1}.
```

单位矩阵不改变向量，也不改变有向体积，所以缩放因子为 1。

### 2. 交换任意两行，determinant 变号

若 $`B`$ 由 $`A`$ 交换任意两行得到，则

```math
\boxed{\det(B)=-\det(A)}.
```

这里不要求交换相邻行。若一共交换 $`s`$ 次，则

```math
\det(B)=(-1)^s\det(A).
```

因此：

- 偶数次交换不改变最终符号；
- 奇数次交换使最终结果变号。

### 3. determinant 对每一行分别线性

固定其他行，只改变其中一行时，determinant 对这一行满足线性性。

若某一行乘以 $`c`$，则

```math
\det
\begin{bmatrix}
\vdots\\
cR_i\\
\vdots
\end{bmatrix}
=
c\det
\begin{bmatrix}
\vdots\\
R_i\\
\vdots
\end{bmatrix}.
```

若某一行是 $`u+v`$，则

```math
\det
\begin{bmatrix}
\vdots\\
u+v\\
\vdots
\end{bmatrix}
=
\det
\begin{bmatrix}
\vdots\\
u\\
\vdots
\end{bmatrix}
+
\det
\begin{bmatrix}
\vdots\\
v\\
\vdots
\end{bmatrix}.
```

拆分某一行时，其他所有行必须保持相同。这种“对每一行分别线性”的性质称为 **multilinearity（多重线性）**。对列也有完全相同的结论。

## 三、由基本性质推出行操作规则

### 1. 两行相同，则 determinant 为零

假设 $`A`$ 有两行完全相同。交换这两行后矩阵没有变化，因此 determinant 也应保持不变；但交换性质又要求它变号：

```math
\det(A)=-\det(A).
```

所以

```math
2\det(A)=0,
```

即

```math
\boxed{\det(A)=0}.
```

### 2. 两行成倍数或有一行全零，则 determinant 为零

若

```math
R_j=cR_i,
```

则从第 $`j`$ 行提出 $`c`$ 后，矩阵出现两行相同，所以

```math
\boxed{\det(A)=0}.
```

零行是这个结论的特殊情况。若某一行全为零，也可以直接从这一行提出系数 0，因此 determinant 为 0。

这已经把 determinant 与线性相关联系起来：

```math
\boxed{
\text{行或列线性相关}
\quad\Longrightarrow\quad
\det(A)=0
}.
```

对于方阵，反过来也成立，后面会得到完整的等价关系。

### 3. 一行加上另一行的倍数，determinant 不变

考虑行替换

```math
R_j\leftarrow R_j+cR_i,
\qquad i\neq j.
```

利用对第 $`j`$ 行的线性性：

```math
\det
\begin{bmatrix}
\vdots\\
R_i\\
\vdots\\
R_j+cR_i\\
\vdots
\end{bmatrix}
=
\det
\begin{bmatrix}
\vdots\\
R_i\\
\vdots\\
R_j\\
\vdots
\end{bmatrix}
+c
\det
\begin{bmatrix}
\vdots\\
R_i\\
\vdots\\
R_i\\
\vdots
\end{bmatrix}.
```

第二个 determinant 有两行相同，所以等于 0。于是

```math
\boxed{
R_j\leftarrow R_j+cR_i
\quad\Longrightarrow\quad
\det\text{ 不变}
}.
```

### 4. 三种初等行操作的影响

| 行操作 | determinant 的变化 |
| --- | --- |
| $`R_i\leftrightarrow R_j`$ | 乘 $`-1`$ |
| $`R_i\leftarrow cR_i`$ | 乘 $`c`$ |
| $`R_j\leftarrow R_j+cR_i`$，$`i\neq j`$ | 不变 |

这张表是用消元计算 determinant 时必须全程追踪的规则。

## 四、三角矩阵与高斯消元

### 1. 三角矩阵的 determinant

若 $`U`$ 是上三角矩阵：

```math
U=
\begin{bmatrix}
d_1&*&\cdots&*\\
0&d_2&\cdots&*\\
\vdots&\ddots&\ddots&\vdots\\
0&\cdots&0&d_n
\end{bmatrix},
```

则

```math
\boxed{
\det(U)=d_1d_2\cdots d_n
}.
```

下三角矩阵也完全相同：determinant 等于主对角线元素的乘积。

这里的 $`d_i`$ 是三角矩阵的对角元素。若 $`U`$ 是消元得到的阶梯形矩阵且 $`A`$ 可逆，它们也就是非零 pivots。

### 2. 用消元计算 determinant

考虑

```math
A=
\begin{bmatrix}
2&1\\
4&3
\end{bmatrix}.
```

执行

```math
R_2\leftarrow R_2-2R_1
```

得到

```math
U=
\begin{bmatrix}
2&1\\
0&1
\end{bmatrix}.
```

这一步只是把另一行的倍数加到第 2 行，所以 determinant 不变：

```math
\det(A)=\det(U)=2\cdot1=2.
```

若消元过程中交换了 $`s`$ 次行，并且没有进行行缩放，则

```math
\boxed{
\det(A)=(-1)^s\prod_{i=1}^n d_i
}.
```

其中 $`d_i`$ 是最终上三角矩阵的对角元素。

### 3. 为什么不能无条件说“determinant 等于 pivot 乘积”

必须检查消元过程中做过什么：

- 行替换不改变 determinant；
- 每次换行要补一个负号；
- 若把某行乘以 $`c`$，必须把这个因子计入；
- 若把 pivot 行除以 pivot，使 pivot 变成 1，也必须把这个缩放还原。

因此，更准确的说法是：

> **当消元只使用行替换和换行时，$`\det(A)`$ 等于换行符号乘以上三角矩阵的对角线乘积。**

## 五、singular 与 invertible 的分类

对 $`n\times n`$ 方阵，只有两种可能：

```text
A singular：不可逆
A nonsingular / invertible：可逆
```

### 1. singular 情况

若 $`A`$ singular，则

```math
\mathrm{rank}(A)<n.
```

消元后不足 $`n`$ 个 pivots，所以阶梯形中一定出现零行：

```math
U=
\begin{bmatrix}
*&*&*\\
0&*&*\\
0&0&0
\end{bmatrix}.
```

三角矩阵至少有一个对角元素为 0，因此

```math
\det(U)=0.
```

换行只能改变符号，行替换不改变 determinant；无论怎样都不能把 0 变成非零数。所以

```math
\boxed{
A\text{ singular}
\quad\Longrightarrow\quad
\det(A)=0
}.
```

### 2. invertible 情况

若 $`A`$ 可逆，则

```math
\mathrm{rank}(A)=n.
```

消元时会得到 $`n`$ 个非零 pivots：

```math
U=
\begin{bmatrix}
d_1&*&*\\
0&d_2&*\\
0&0&d_3
\end{bmatrix},
\qquad
d_i\neq0.
```

所以

```math
\det(A)
=(-1)^s d_1d_2\cdots d_n
\neq0.
```

因此

```math
\boxed{
A\text{ invertible}
\quad\Longrightarrow\quad
\det(A)\neq0
}.
```

singular 与 invertible 已经覆盖所有方阵：前者的 determinant 为 0，后者的 determinant 非零。因此

```math
\boxed{
A\text{ singular}
\quad\Longleftrightarrow\quad
\det(A)=0
}.
```

### 3. 与前面知识的完整等价链

对 $`n\times n`$ 方阵 $`A`$：

```math
\boxed{
\begin{aligned}
A\text{ invertible}
&\Longleftrightarrow \det(A)\neq0\\
&\Longleftrightarrow \mathrm{rank}(A)=n\\
&\Longleftrightarrow A\text{ 有 }n\text{ 个 pivots}\\
&\Longleftrightarrow N(A)=\{\mathbf0\}\\
&\Longleftrightarrow A\text{ 的列线性无关}\\
&\Longleftrightarrow \mathrm{Col}(A)=\mathbb{R}^n\\
&\Longleftrightarrow A\mathbf{x}=\mathbf{b}
\text{ 对每个 }\mathbf{b}\text{ 都有唯一解}.
\end{aligned}
}
```

singular 情况就是这些条件的反面：

```math
\boxed{
A\text{ singular}
\Longleftrightarrow
\det(A)=0
\Longleftrightarrow
\mathrm{rank}(A)<n
\Longleftrightarrow
N(A)\neq\{\mathbf0\}
}.
```

需要注意：当 $`A`$ singular 时，$`A\mathbf{x}=\mathbf{b}`$ 可能无解，也可能有无穷多解；它不可能对每个 $`\mathbf{b}`$ 都有唯一解。

> 老师的 case split 可以概括为：singular 时消元会出现零行；invertible 时会得到 $`n`$ 个非零 pivots，其乘积给出非零 determinant。

## 六、LU 分解来自高斯消元

### 1. $`A=LU`$ 不是 determinant 的定义

```math
\boxed{A=LU}
```

是矩阵分解，不是 determinant 的公理。它把高斯消元

```math
A\longrightarrow U
```

中使用的消元倍数记录到下三角矩阵 $`L`$ 中：

- $`U`$：消元后的上三角矩阵；
- $`L`$：记录 elimination multipliers 的单位下三角矩阵。

### 2. 从 elimination matrix 推出 $`A=LU`$

仍取

```math
A=
\begin{bmatrix}
2&1\\
4&3
\end{bmatrix}.
```

行操作

```math
R_2\leftarrow R_2-2R_1
```

可以写成左乘 elimination matrix：

```math
E=
\begin{bmatrix}
1&0\\
-2&1
\end{bmatrix},
\qquad
EA=U.
```

这里

```math
U=
\begin{bmatrix}
2&1\\
0&1
\end{bmatrix}.
```

由 $`EA=U`$ 得

```math
A=E^{-1}U.
```

定义

```math
L=E^{-1}
=
\begin{bmatrix}
1&0\\
2&1
\end{bmatrix},
```

于是

```math
\boxed{A=LU}.
```

直接验证：

```math
\begin{bmatrix}
1&0\\
2&1
\end{bmatrix}
\begin{bmatrix}
2&1\\
0&1
\end{bmatrix}
=
\begin{bmatrix}
2&1\\
4&3
\end{bmatrix}.
```

$`L`$ 左下角的 2，正是消元时使用的倍数。

### 3. 多步消元

若连续消元满足

```math
E_k\cdots E_2E_1A=U,
```

则

```math
A
=E_1^{-1}E_2^{-1}\cdots E_k^{-1}U.
```

把前面的逆矩阵乘积记为

```math
L=E_1^{-1}E_2^{-1}\cdots E_k^{-1},
```

便得到

```math
\boxed{A=LU}.
```

这说明：

```math
\boxed{
LU\text{ 分解是 Gaussian elimination 的矩阵化记录}
}.
```

### 4. 需要换行时：$`PA=LU`$

并不是每个矩阵都能在完全不换行的情况下直接得到普通的 $`A=LU`$。例如

```math
A=
\begin{bmatrix}
0&1\\
1&2
\end{bmatrix}
```

的第一个候选 pivot 是 0，必须先交换两行。令

```math
P=
\begin{bmatrix}
0&1\\
1&0
\end{bmatrix},
```

则更一般的分解写成

```math
\boxed{PA=LU}.
```

$`P`$ 是 permutation matrix（置换矩阵），记录消元过程中的换行。实际数值计算通常使用这种带 pivoting 的形式。

### 5. LU 怎样帮助计算 determinant

在不换行的 $`A=LU`$ 中，$`L`$ 通常是单位下三角矩阵，所以

```math
\det(L)=1.
```

因此

```math
\det(A)
=\det(L)\det(U)
=\det(U)
=\prod_{i=1}^n u_{ii}.
```

若

```math
PA=LU
```

并且 $`P`$ 记录了 $`s`$ 次换行，则

```math
\det(P)=(-1)^s.
```

由

```math
\det(P)\det(A)=\det(L)\det(U)
```

得到

```math
\boxed{
\det(A)=(-1)^s\prod_{i=1}^n u_{ii}
}.
```

所以逻辑顺序是：先对矩阵做 LU 分解，再对等式两边取 determinant；不是“把 determinant 分解成 LU”。

## 七、determinant 的代数性质

### 1. 乘积的 determinant

对同阶方阵 $`A,B`$：

```math
\boxed{
\det(AB)=\det(A)\det(B)
}.
```

这表示连续做两个线性变换时，总体积缩放因子等于两个缩放因子的乘积。

矩阵乘法一般不满足 $`AB=BA`$，但 determinant 是数，因此

```math
\det(AB)
=\det(A)\det(B)
=\det(B)\det(A)
=\det(BA).
```

所以即使 $`AB\neq BA`$，仍可能有

```math
\det(AB)=\det(BA).
```

### 2. 幂与逆矩阵

由乘法性质：

```math
\det(A^2)
=\det(A)\det(A)
=\bigl(\det(A)\bigr)^2.
```

更一般地，正整数 $`k`$ 满足

```math
\det(A^k)=\bigl(\det(A)\bigr)^k.
```

若 $`A`$ 可逆，则

```math
1
=\det(I)
=\det(AA^{-1})
=\det(A)\det(A^{-1}),
```

因此

```math
\boxed{
\det(A^{-1})=\frac{1}{\det(A)}
}.
```

### 3. 整个矩阵乘常数

若 $`A`$ 是 $`n\times n`$ 方阵，则 $`cA`$ 把 $`A`$ 的每一行都乘以 $`c`$。一共有 $`n`$ 行，所以

```math
\boxed{
\det(cA)=c^n\det(A)
}.
```

例如，对于 $`3\times3`$ 矩阵：

```math
\det(2A)=2^3\det(A)=8\det(A).
```

必须区分：

- 只把一行乘 $`c`$：determinant 乘 $`c`$；
- 把整个 $`n\times n`$ 矩阵乘 $`c`$：determinant 乘 $`c^n`$。

### 4. 转置不改变 determinant

```math
\boxed{
\det(A^T)=\det(A)
}.
```

这条性质说明 determinant 对行与列是对称的。前面对行得到的规则，也都可以改写成对应的列规则。

可以按 singular 与 invertible 分两种情况理解。

#### 情况一：$`A`$ singular

因为

```math
\mathrm{rank}(A^T)=\mathrm{rank}(A)<n,
```

所以 $`A^T`$ 也 singular。因此

```math
\det(A)=\det(A^T)=0.
```

#### 情况二：$`A`$ invertible

用带换行的 LU 分解写成

```math
PA=LU.
```

由于 $`P^{-1}=P^T`$，有

```math
A=P^TLU.
```

转置后：

```math
A^T=U^TL^TP.
```

上三角矩阵转置后变成下三角矩阵，但主对角线不变，所以

```math
\det(U^T)=\det(U),
\qquad
\det(L^T)=\det(L).
```

置换矩阵 $`P`$ 与 $`P^T`$ 记录相同的交换奇偶性，因此

```math
\det(P^T)=\det(P).
```

于是

```math
\begin{aligned}
\det(A^T)
&=\det(U^T)\det(L^T)\det(P)\\
&=\det(U)\det(L)\det(P)\\
&=\det(P^T)\det(L)\det(U)\\
&=\det(A).
\end{aligned}
```

若消元不需要换行，$`P=I`$，上面的证明就简化为老师板书中的

```math
A=LU,
\qquad
A^T=U^TL^T.
```

## 八、常见混淆与修正

### 1. “行列式可以分解成 LU”

不准确。正确顺序是：

```text
矩阵 A 经过消元得到矩阵分解 A=LU 或 PA=LU
    -> 再对矩阵等式两边取 determinant
    -> 得到 det(A) 与 det(L)、det(U) 的关系
```

### 2. “任何矩阵都有 $`A=LU`$”

不准确。若消元不需要换行，可以直接写 $`A=LU`$；需要换行时，一般写成 $`PA=LU`$。

### 3. “交换远距离的两行要按相邻交换次数算”

如果直接执行一次 $`R_i\leftrightarrow R_j`$，这就是一次行交换，determinant 乘 $`-1`$。若人为把它拆成若干次相邻交换，则必须计算整个置换的奇偶性；最终符号与直接交换一致。

### 4. “只要做消元，determinant 就是对角线乘积”

必须追踪换行与行缩放。只有行替换不改变 determinant。

### 5. “$`\det(2A)=2\det(A)`$”

对 $`n\times n`$ 方阵应为

```math
\det(2A)=2^n\det(A).
```

因为整个矩阵的每一行都乘了 2。

## 九、本章总结

### 1. 最核心的计算链

```text
A
    -> 高斯消元
U
    -> 追踪换行次数与行缩放
det(A)=符号与缩放修正 × U 的对角线乘积
```

若消元只使用行替换和 $`s`$ 次换行：

```math
\boxed{
\det(A)=(-1)^s\prod_{i=1}^n u_{ii}
}.
```

### 2. 最核心的结构链

```math
\boxed{
\det(A)=0
\Longleftrightarrow
A\text{ singular}
\Longleftrightarrow
\mathrm{rank}(A)<n
\Longleftrightarrow
\text{少于 }n\text{ 个 pivots}
\Longleftrightarrow
N(A)\neq\{\mathbf0\}
\Longleftrightarrow
\text{columns dependent}
}.
```

反过来：

```math
\boxed{
\det(A)\neq0
\Longleftrightarrow
A\text{ invertible}
\Longleftrightarrow
\mathrm{rank}(A)=n
\Longleftrightarrow
N(A)=\{\mathbf0\}
}.
```

### 3. 最值得记住的公式

```math
\boxed{det(I)=1}
```

```math
\boxed{det(U)=\prod_{i=1}^n u_{ii}}
```

```math
\boxed{det(AB)=\det(A)\det(B)}
```

```math
\boxed{det(A^T)=\det(A)}
```

```math
\boxed{det(cA)=c^n\det(A)}
```

```math
\boxed{det(A^{-1})=\frac{1}{\det(A)}}
\qquad(A\text{ invertible}).
```

### 4. 本章终点

determinant 不是与前面内容割裂的新计算技巧。它把一个方阵的可逆性、满秩、pivot 数量、零空间、列独立和体积缩放统一到一个标量上：

```text
det(A) 是否为 0
    -> 判断 A 是否丢失维度
    -> 判断 A 是否 singular
    -> 判断列是否独立
    -> 判断 Ax=b 是否对每个 b 都有唯一解
```

而 LU 分解则解释了为什么 determinant 能通过消元高效计算：$`L`$ 记录消元过程，$`U`$ 暴露 pivots 与对角线乘积。
