# Review 1 习题

本文件用于整理 Review 1 的习题。每道题按照“题目条件、核心思路、完整推导、验证与方法总结”的顺序记录。

## 第 1 题：由完整解反推矩阵

### 1. 题目条件

已知 $`A`$ 是一个 $`3\times3`$ 矩阵，方程

```math
A\mathbf{x}
=
\begin{bmatrix}
2\\4\\2
\end{bmatrix}
```

的全部解为

```math
\mathbf{x}
=
\begin{bmatrix}
2\\0\\0
\end{bmatrix}
+c
\begin{bmatrix}
1\\1\\0
\end{bmatrix}
+d
\begin{bmatrix}
0\\0\\1
\end{bmatrix},
\qquad c,d\in\mathbb{R}.
```

要求根据完整解确定矩阵 $`A`$，并说明它的零空间、列关系与秩。

### 2. 先识别完整解的结构

非齐次线性方程的全部解具有统一形式

```math
\mathbf{x}=\mathbf{x}_p+\mathbf{x}_n,
\qquad
\mathbf{x}_n\in N(A).
```

题目给出的第一个向量是一个特解：

```math
\mathbf{x}_p
=
\begin{bmatrix}
2\\0\\0
\end{bmatrix}.
```

后面带有任意参数 $`c,d`$ 的部分是零空间中的任意向量：

```math
\mathbf{x}_n
=c\mathbf{v}_1+d\mathbf{v}_2,
```

其中

```math
\mathbf{v}_1
=
\begin{bmatrix}
1\\1\\0
\end{bmatrix},
\qquad
\mathbf{v}_2
=
\begin{bmatrix}
0\\0\\1
\end{bmatrix}.
```

因此

```math
N(A)
=\mathrm{span}
\left\{
\begin{bmatrix}1\\1\\0\end{bmatrix},
\begin{bmatrix}0\\0\\1\end{bmatrix}
\right\}.
```

这两个向量线性无关，所以它们构成 $`N(A)`$ 的一组基：

```math
\dim N(A)=2.
```

> **读完整解时，第一个任务不是立刻计算矩阵，而是先分清哪一部分负责到达目标输出，哪一部分只在零空间中移动。**

### 3. 为什么加入零空间向量不会改变输出

记

```math
\mathbf{b}
=
\begin{bmatrix}
2\\4\\2
\end{bmatrix}.
```

特解与零空间向量分别满足

```math
A\mathbf{x}_p=\mathbf{b},
\qquad
A\mathbf{v}_1=\mathbf{0},
\qquad
A\mathbf{v}_2=\mathbf{0}.
```

利用线性性，对任意 $`c,d`$ 都有

```math
A(\mathbf{x}_p+c\mathbf{v}_1+d\mathbf{v}_2)
=A\mathbf{x}_p+cA\mathbf{v}_1+dA\mathbf{v}_2
=\mathbf{b}.
```

所以完整解中的两项作用不同：

| 部分 | 作用 |
| --- | --- |
| $`\mathbf{x}_p`$ | 把输出送到目标 $`\mathbf{b}`$ |
| $`c\mathbf{v}_1+d\mathbf{v}_2`$ | 在零空间中自由移动，不改变 $`A\mathbf{x}`$ |

### 4. 按列反推矩阵

把 $`A`$ 的三列记为

```math
A=
\begin{bmatrix}
\mathbf{a}_1&\mathbf{a}_2&\mathbf{a}_3
\end{bmatrix}.
```

矩阵乘向量就是用向量的三个分量作为系数，对 $`A`$ 的三列作线性组合。下面分别使用特解和两个零空间基向量。

#### 第一步：由特解确定第一列

因为

```math
A
\begin{bmatrix}
2\\0\\0
\end{bmatrix}
=
\begin{bmatrix}
2\\4\\2
\end{bmatrix},
```

按列展开可得

```math
2\mathbf{a}_1
=
\begin{bmatrix}
2\\4\\2
\end{bmatrix}.
```

所以

```math
\mathbf{a}_1
=
\begin{bmatrix}
1\\2\\1
\end{bmatrix}.
```

#### 第二步：由第一个零空间方向确定第二列

因为 $`\mathbf{v}_1\in N(A)`$，所以

```math
A
\begin{bmatrix}
1\\1\\0
\end{bmatrix}
=\mathbf{0}.
```

按列展开：

```math
\mathbf{a}_1+\mathbf{a}_2=\mathbf{0}.
```

因此

```math
\mathbf{a}_2
=-\mathbf{a}_1
=
\begin{bmatrix}
-1\\-2\\-1
\end{bmatrix}.
```

#### 第三步：由第二个零空间方向确定第三列

因为 $`\mathbf{v}_2\in N(A)`$，所以

```math
A
\begin{bmatrix}
0\\0\\1
\end{bmatrix}
=\mathbf{0}.
```

按列展开：

```math
\mathbf{a}_3=\mathbf{0}.
```

因此

```math
\mathbf{a}_3
=
\begin{bmatrix}
0\\0\\0
\end{bmatrix}.
```

三列合在一起，得到

```math
\boxed{
A=
\begin{bmatrix}
1&-1&0\\
2&-2&0\\
1&-1&0
\end{bmatrix}
}.
```

### 5. 为什么这个矩阵能够被唯一确定

题目给出的三个输入向量

```math
\mathbf{x}_p,
\qquad
\mathbf{v}_1,
\qquad
\mathbf{v}_2
```

线性无关，因此构成 $`\mathbb{R}^3`$ 的一组基。题目同时告诉了 $`A`$ 对这三个基向量分别产生什么输出：

```math
A\mathbf{x}_p=\mathbf{b},
\qquad
A\mathbf{v}_1=\mathbf{0},
\qquad
A\mathbf{v}_2=\mathbf{0}.
```

线性变换在一组基上的输出一旦确定，它对整个空间的作用也就确定了，所以这里得到的 $`A`$ 是唯一的。

### 6. 验证全部解

先把完整解合并：

```math
\mathbf{x}
=
\begin{bmatrix}
2+c\\c\\d
\end{bmatrix}.
```

代入求出的矩阵：

```math
A\mathbf{x}
=
\begin{bmatrix}
1&-1&0\\
2&-2&0\\
1&-1&0
\end{bmatrix}
\begin{bmatrix}
2+c\\c\\d
\end{bmatrix}.
```

逐行计算得到

```math
A\mathbf{x}
=
\begin{bmatrix}
(2+c)-c\\
2(2+c)-2c\\
(2+c)-c
\end{bmatrix}
=
\begin{bmatrix}
2\\4\\2
\end{bmatrix}.
```

结果与 $`c,d`$ 无关，因此题目给出的所有向量确实都是原方程的解。

### 7. 从零空间判断秩

$`A`$ 有三个未知数，所以列数为 $`n=3`$。由秩-零度关系

```math
\dim N(A)=n-r
```

以及 $`\dim N(A)=2`$，得到

```math
2=3-r,
\qquad
r=1.
```

所以

```math
\mathrm{rank}(A)=1.
```

这与矩阵的列关系完全一致：

```math
\mathbf{a}_2=-\mathbf{a}_1,
\qquad
\mathbf{a}_3=\mathbf{0}.
```

三列中只有 $`\mathbf{a}_1`$ 提供一个独立方向，因此列空间的一组基可以取为

```math
\left\{
\begin{bmatrix}1\\2\\1\end{bmatrix}
\right\},
```

列空间维数与矩阵的秩都是 $`1`$。

### 8. 解题逻辑链

```text
读出完整解 x = x_p + c v_1 + d v_2
    -> x_p 是特解，A x_p = b
    -> v_1、v_2 是零空间方向，A v_1 = A v_2 = 0
    -> 按列展开三条矩阵方程
    -> 得到 a_1、a_2、a_3
    -> 组合出矩阵 A
    -> 用 dim N(A) = n-r 得到 rank(A)=1
    -> 用 a_2=-a_1、a_3=0 验证只有一个独立列方向
```

### 9. 最终答案

```math
\boxed{
A=
\begin{bmatrix}
1&-1&0\\
2&-2&0\\
1&-1&0
\end{bmatrix}
}
```

```math
\boxed{
N(A)
=\mathrm{span}
\left\{
\begin{bmatrix}1\\1\\0\end{bmatrix},
\begin{bmatrix}0\\0\\1\end{bmatrix}
\right\},
\qquad
\dim N(A)=2
}
```

```math
\boxed{
\mathrm{rank}(A)=1,
\qquad
\mathbf{a}_2=-\mathbf{a}_1,
\qquad
\mathbf{a}_3=\mathbf{0}
}
```

> **这道题的核心方法：完整解不仅描述所有 $`\mathbf{x}`$，还同时告诉了线性变换 $`A`$ 在特解方向和全部零空间方向上的作用。把这些信息按列展开，就能反推出矩阵，并进一步读出列关系与秩。**

### 10. 对哪些右端向量方程有解

求出 $`A`$ 后，它的三列满足

```math
\mathbf{a}_1
=
\begin{bmatrix}1\\2\\1\end{bmatrix},
\qquad
\mathbf{a}_2=-\mathbf{a}_1,
\qquad
\mathbf{a}_3=\mathbf{0}.
```

因此所有列只有一个独立方向，列空间为

```math
\mathrm{Col}(A)
=\mathrm{span}
\left\{
\begin{bmatrix}1\\2\\1\end{bmatrix}
\right\}.
```

方程是否有解由右端向量是否属于列空间决定：

```math
A\mathbf{x}=\mathbf{b}\text{ 有解}
\quad\Longleftrightarrow\quad
\mathbf{b}\in\mathrm{Col}(A).
```

所以本题中可解的右端向量必须具有形式

```math
\boxed{
\mathbf{b}
=k
\begin{bmatrix}1\\2\\1\end{bmatrix},
\qquad k\in\mathbb{R}
}.
```

直接计算也能看出同一个结论。对任意

```math
\mathbf{x}
=
\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix},
```

都有

```math
A\mathbf{x}
=
\begin{bmatrix}
x_1-x_2\\
2x_1-2x_2\\
x_1-x_2
\end{bmatrix}.
```

令 $`k=x_1-x_2`$，便得到

```math
A\mathbf{x}
=k
\begin{bmatrix}1\\2\\1\end{bmatrix}.
```

无论怎样选择 $`\mathbf{x}`$，输出都只能落在 $`(1,2,1)^T`$ 张成的直线上。例如，$`(1,2,1)^T`$、$`(3,6,3)^T`$ 和 $`(-2,-4,-2)^T`$ 都可以产生，而

```math
\begin{bmatrix}1\\0\\0\end{bmatrix}
\notin\mathrm{Col}(A),
```

所以以它为右端向量时方程无解。

如果

```math
\mathbf{b}
=k
\begin{bmatrix}1\\2\\1\end{bmatrix},
```

可以取一个特解

```math
\mathbf{x}_p
=
\begin{bmatrix}k\\0\\0\end{bmatrix}.
```

于是对应的全部解为

```math
\mathbf{x}
=
\begin{bmatrix}k\\0\\0\end{bmatrix}
+c\begin{bmatrix}1\\1\\0\end{bmatrix}
+d\begin{bmatrix}0\\0\\1\end{bmatrix}.
```

原题的右端向量是 $`2(1,2,1)^T`$，所以取 $`k=2`$ 就回到了原来的完整解。

### 11. 从本例推广到存在性与唯一性

设

```math
A\in\mathbb{R}^{m\times n},
\qquad
r=\mathrm{rank}(A).
```

**满行秩**是指 $`r=m`$。此时

```math
\dim\bigl(\mathrm{Col}(A)\bigr)=m.
```

因为列空间本来就是 $`\mathbb{R}^m`$ 的子空间，一个 $`m`$ 维子空间只能是整个 $`\mathbb{R}^m`$：

```math
\mathrm{Col}(A)=\mathbb{R}^m.
```

所以每一个 $`\mathbf{b}\in\mathbb{R}^m`$ 都能到达。

> **$`r=m`$ 控制存在性：对每个右端向量，$`A\mathbf{x}=\mathbf{b}`$ 都有解。**

**满列秩**是指 $`r=n`$。由秩-零度关系，

```math
\dim N(A)=n-r=0,
\qquad
N(A)=\{\mathbf{0}\}.
```

如果同一个右端向量有两个解 $`\mathbf{x}^{(1)}`$ 和 $`\mathbf{x}^{(2)}`$，那么

```math
A(\mathbf{x}^{(1)}-\mathbf{x}^{(2)})=\mathbf{0}.
```

零空间只有零向量，所以两解之差只能为零，两解必须相同。

> **$`r=n`$ 控制唯一性：方程一旦有解，这个解就是唯一的。**

两条结论不能混用：满行秩保证“至少一个解”，满列秩保证“至多一个解”。只有二者同时成立，才会得到“恰好一个解”。

| 秩条件 | 列空间 | 零空间 | 对 $`A\mathbf{x}=\mathbf{b}`$ 的影响 |
| --- | --- | --- | --- |
| $`r=m`$ | $`\mathrm{Col}(A)=\mathbb{R}^m`$ | 不一定只有零向量 | 每个 $`\mathbf{b}`$ 至少有一个解 |
| $`r=n`$ | 不一定等于 $`\mathbb{R}^m`$ | $`N(A)=\{\mathbf{0}\}`$ | 有解时至多一个解 |
| $`r=m=n`$ | $`\mathrm{Col}(A)=\mathbb{R}^n`$ | $`N(A)=\{\mathbf{0}\}`$ | 每个 $`\mathbf{b}`$ 恰好有一个解 |

### 12. 方阵中的零空间与可逆性

若 $`A`$ 是 $`n\times n`$ 方阵，并且

```math
N(A)=\{\mathbf{0}\},
```

那么 $`\dim N(A)=0`$。由

```math
\dim N(A)=n-r
```

得到 $`r=n`$。方阵又满足 $`m=n`$，所以

```math
r=m=n.
```

转置不改变秩：

```math
\mathrm{rank}(A^T)=\mathrm{rank}(A)=n.
```

$`A^T`$ 仍然是 $`n\times n`$ 方阵，因此

```math
\dim N(A^T)
=n-\mathrm{rank}(A^T)
=0,
```

从而

```math
N(A^T)=\{\mathbf{0}\}.
```

这条推理可以压缩为

```math
N(A)=\{\mathbf{0}\}
\quad\Longrightarrow\quad
r=n=m
\quad\Longrightarrow\quad
N(A^T)=\{\mathbf{0}\}.
```

对于方阵，下面这些说法描述的是同一个现象：

```text
列线性无关
    <-> N(A)={0}
    <-> rank(A)=n
    <-> Col(A)=R^n
    <-> Ax=b 对每个 b 都有唯一解
    <-> A 可逆
    <-> N(A^T)={0}
```

这里“$`A`$ 是方阵”是不可缺少的条件。例如，

```math
A=
\begin{bmatrix}
1&0\\
0&1\\
0&0
\end{bmatrix}
\in\mathbb{R}^{3\times2}
```

具有满列秩，所以 $`N(A)=\{\mathbf{0}\}`$；但是

```math
A^T=
\begin{bmatrix}
1&0&0\\
0&1&0
\end{bmatrix}
```

满足

```math
N(A^T)
=\mathrm{span}
\left\{
\begin{bmatrix}0\\0\\1\end{bmatrix}
\right\}
\neq\{\mathbf{0}\}.
```

矩形矩阵中，$`r=n`$ 只保证满列秩和唯一性，不一定同时有 $`r=m`$；方阵中 $`m=n`$，所以满列秩会自动成为满行秩，同时保证存在性与唯一性。

### 13. 第 1 题的完整知识链

```text
完整解 x=x_p+x_n
    -> 从 x_p 得到一个可达右端向量
    -> 从 N(A) 的基得到矩阵列之间的关系
    -> 反推出 A，并由 nullity 得到 rank(A)=1
    -> Col(A) 是由 (1,2,1)^T 张成的一条直线
    -> Ax=b 有解当且仅当 b 在这条直线上
    -> 推广：r=m 控制存在性，r=n 控制唯一性
    -> 方阵中二者会同时成立，等价于 A 可逆
```

## 第 2 题：利用可逆左乘因子求零空间

### 1. 题目条件

已知

```math
B=CD,
```

其中

```math
C=
\begin{bmatrix}
1&1&0\\
0&1&0\\
1&0&1
\end{bmatrix},
\qquad
D=
\begin{bmatrix}
1&0&-1&2\\
0&1&1&-1\\
0&0&0&0
\end{bmatrix}.
```

要求求出 $`N(B)`$ 的一组基和维数。

矩阵尺寸为

```math
C\in\mathbb{R}^{3\times3},
\qquad
D\in\mathbb{R}^{3\times4},
\qquad
B\in\mathbb{R}^{3\times4}.
```

所以

```math
N(B)\subseteq\mathbb{R}^4.
```

这只表示零空间中的向量具有四个分量，不表示 $`N(B)=\mathbb{R}^4`$。

### 2. 为什么可以直接研究 $`D`$

先计算

```math
\det(C)=1\neq0,
```

所以 $`C`$ 可逆。对任意 $`\mathbf{x}\in\mathbb{R}^4`$，

```math
B\mathbf{x}=\mathbf{0}
\quad\Longleftrightarrow\quad
CD\mathbf{x}=\mathbf{0}.
```

由于 $`C`$ 可逆，可以在等式左边乘 $`C^{-1}`$：

```math
CD\mathbf{x}=\mathbf{0}
\quad\Longleftrightarrow\quad
D\mathbf{x}=\mathbf{0}.
```

因此

```math
\boxed{N(B)=N(CD)=N(D)}.
```

直观上，输入先经过 $`D`$，再经过 $`C`$：

```math
\mathbf{x}
\xrightarrow{\ D\ }
D\mathbf{x}
\xrightarrow{\ C\ }
CD\mathbf{x}.
```

可逆矩阵 $`C`$ 不会把一个非零向量压成零向量，所以 $`CD\mathbf{x}=\mathbf{0}`$ 与 $`D\mathbf{x}=\mathbf{0}`$ 具有完全相同的解。

> **左乘可逆矩阵会改变输出的表示，但不会改变哪些输入被送到零向量。**

若 $`C`$ 不可逆，只能保证 $`N(D)\subseteq N(CD)`$，因为某些非零的 $`D\mathbf{x}`$ 还可能被 $`C`$ 进一步压成零；因此“$`C`$ 可逆”是等号成立的关键条件。

### 3. 从 $`D`$ 求零空间

设

```math
\mathbf{x}
=
\begin{bmatrix}
x_1\\x_2\\x_3\\x_4
\end{bmatrix}.
```

方程 $`D\mathbf{x}=\mathbf{0}`$ 给出

```math
x_1-x_3+2x_4=0,
\qquad
x_2+x_3-x_4=0.
```

第 $`1,2`$ 列是主元列，所以 $`x_3,x_4`$ 是自由变量。令

```math
x_3=s,
\qquad
x_4=t.
```

主元变量由它们决定：

```math
x_1=s-2t,
\qquad
x_2=-s+t.
```

因此

```math
\mathbf{x}
=
\begin{bmatrix}
s-2t\\-s+t\\s\\t
\end{bmatrix}
=s
\begin{bmatrix}
1\\-1\\1\\0
\end{bmatrix}
+t
\begin{bmatrix}
-2\\1\\0\\1
\end{bmatrix}.
```

由 $`N(B)=N(D)`$，得到

```math
\boxed{
N(B)
=\mathrm{span}
\left\{
\begin{bmatrix}1\\-1\\1\\0\end{bmatrix},
\begin{bmatrix}-2\\1\\0\\1\end{bmatrix}
\right\}
}.
```

这两个特殊解线性无关，因此构成 $`N(B)`$ 的一组基。

### 4. 维数与秩

$`D`$ 有四列和两个主元，所以

```math
\mathrm{rank}(D)=2,
\qquad
\dim N(D)=4-2=2.
```

左乘可逆矩阵也不改变秩，因此

```math
\mathrm{rank}(B)=\mathrm{rank}(D)=2.
```

结合 $`N(B)=N(D)`$，得到

```math
\boxed{\dim N(B)=2}.
```

### 5. 直接验证

先计算

```math
B=CD
=
\begin{bmatrix}
1&1&0&1\\
0&1&1&-1\\
1&0&-1&2
\end{bmatrix}.
```

令

```math
\mathbf{v}_1
=
\begin{bmatrix}1\\-1\\1\\0\end{bmatrix},
\qquad
\mathbf{v}_2
=
\begin{bmatrix}-2\\1\\0\\1\end{bmatrix}.
```

直接相乘可得

```math
B\mathbf{v}_1=\mathbf{0},
\qquad
B\mathbf{v}_2=\mathbf{0}.
```

所以两个向量确实都属于 $`N(B)`$；又因为零空间维数为 $`2`$，它们正好构成一组基。

### 6. 与高斯消元的关系

每一步行操作都可以写成左乘一个可逆的消元矩阵。若一系列行操作把 $`A`$ 化为 $`R`$，就存在可逆矩阵 $`E`$ 使

```math
R=EA.
```

于是

```math
R\mathbf{x}=\mathbf{0}
\quad\Longleftrightarrow\quad
EA\mathbf{x}=\mathbf{0}
\quad\Longleftrightarrow\quad
A\mathbf{x}=\mathbf{0},
```

所以

```math
N(R)=N(A).
```

这正是“可以从行阶梯形或 RREF 求原矩阵零空间”的根本原因。

### 7. 第 2 题的解题逻辑链

```text
观察 B=CD
    -> 检查 C 是否可逆
    -> C 可逆，所以 CDx=0 <-> Dx=0
    -> N(B)=N(D)
    -> 从简单的 D 识别主元变量与自由变量
    -> 写出两个特殊解
    -> 得到 N(B) 的基与 dim N(B)=2
    -> 认出这就是“行操作不改变零空间”的一般原理
```

## Review 1 当前总结

这两道题从两个方向使用了零空间：

1. 第 1 题从“特解 + 零空间”反推矩阵，并由列空间判断哪些右端向量可解。
2. 第 2 题利用可逆左乘不改变零空间，把复杂矩阵的齐次方程转化为更简单矩阵的齐次方程。

它们共同使用的核心关系是

```math
\dim N(A)=n-r,
\qquad
\dim\bigl(\mathrm{Col}(A)\bigr)=r,
```

```math
A\mathbf{x}=\mathbf{b}\text{ 有解}
\quad\Longleftrightarrow\quad
\mathbf{b}\in\mathrm{Col}(A),
```

```math
C\text{ 可逆}
\quad\Longrightarrow\quad
N(CD)=N(D).
```

> **Review 1 的主线不是孤立地计算矩阵，而是从解的结构识别零空间，从零空间识别自由度与列关系，再用列空间、秩和可逆性判断方程的存在性与唯一性。**
