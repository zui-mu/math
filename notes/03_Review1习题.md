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

## 第 3 题：交换两行对四个基本子空间的影响

### 1. 题目条件

交换矩阵 $`A`$ 的两行，相当于左乘一个置换矩阵 $`P`$：

```math
B=PA.
```

置换矩阵是可逆矩阵，并且

```math
P^{-1}=P^T.
```

需要判断交换行后，$`A`$ 的四个基本子空间哪些保持不变，哪些通常会改变。

结论如下：

| 空间 | 交换行后是否保持为原来的空间 |
| --- | --- |
| $`N(A)`$ | 不变 |
| $`\mathrm{Row}(A)=\mathrm{Col}(A^T)`$ | 不变 |
| $`\mathrm{Col}(A)`$ | 通常改变 |
| $`N(A^T)`$ | 通常改变 |

虽然空间本身有的会改变，但四个空间的维数都不变，因为左乘可逆矩阵不改变秩。

### 2. 零空间不变

```math
B\mathbf{x}=\mathbf{0}
\quad\Longleftrightarrow\quad
PA\mathbf{x}=\mathbf{0}.
```

由于 $`P`$ 可逆，可以左乘 $`P^{-1}`$：

```math
PA\mathbf{x}=\mathbf{0}
\quad\Longleftrightarrow\quad
A\mathbf{x}=\mathbf{0}.
```

所以

```math
\boxed{N(B)=N(A)}.
```

这就是行交换以及一般初等行变换不改变齐次方程解集的原因。

### 3. 行空间不变

若 $`A`$ 的行向量为 $`\mathbf{r}_1,\mathbf{r}_2,\ldots`$，交换两行只是改变这些行的排列顺序：

```math
\mathrm{span}\{\mathbf{r}_1,\mathbf{r}_2,\ldots\}
=
\mathrm{span}\{\mathbf{r}_2,\mathbf{r}_1,\ldots\}.
```

因此

```math
\boxed{\mathrm{Row}(B)=\mathrm{Row}(A)}
```

也就是

```math
\boxed{\mathrm{Col}(B^T)=\mathrm{Col}(A^T)}.
```

更一般地，初等行变换产生的新行是旧行的可逆线性组合，所以行空间保持不变。

### 4. 列空间通常会改变

把 $`A`$ 按列写成

```math
A=\begin{bmatrix}\mathbf{a}_1&\cdots&\mathbf{a}_n\end{bmatrix},
```

那么

```math
B=PA
=\begin{bmatrix}P\mathbf{a}_1&\cdots&P\mathbf{a}_n\end{bmatrix}.
```

因此

```math
\boxed{\mathrm{Col}(B)=P\,\mathrm{Col}(A)}.
```

这表示列空间整体经过 $`P`$ 的坐标置换。它的维数不变，但在输出空间中的位置可能改变。

例如

```math
A=\begin{bmatrix}1&0\\0&0\end{bmatrix},
\qquad
B=\begin{bmatrix}0&0\\1&0\end{bmatrix},
```

其中 $`B`$ 是交换 $`A`$ 的两行得到的矩阵。此时

```math
\mathrm{Col}(A)
=\mathrm{span}\left\{\begin{bmatrix}1\\0\end{bmatrix}\right\},
\qquad
\mathrm{Col}(B)
=\mathrm{span}\left\{\begin{bmatrix}0\\1\end{bmatrix}\right\}.
```

原来的列空间是 $`x`$ 轴方向，交换后变成 $`y`$ 轴方向；但二者维数都为 $`1`$。

因此必须区分：

> **行变换保持秩和列空间的维数，但通常不保持列空间本身。**

### 5. 左零空间也会随输出坐标改变

由 $`B=PA`$ 可得

```math
B^T=A^TP^T.
```

如果 $`\mathbf{y}\in N(B^T)`$，则

```math
A^TP^T\mathbf{y}=\mathbf{0}.
```

令 $`\mathbf{z}=P^T\mathbf{y}`$，就有 $`\mathbf{z}\in N(A^T)`$。又因为 $`\mathbf{y}=P\mathbf{z}`$，所以

```math
\boxed{N(B^T)=P\,N(A^T)}.
```

继续使用上面的例子：

```math
N(A^T)
=\mathrm{span}\left\{\begin{bmatrix}0\\1\end{bmatrix}\right\},
\qquad
N(B^T)
=\mathrm{span}\left\{\begin{bmatrix}1\\0\end{bmatrix}\right\}.
```

所以左零空间也改变了，但维数仍然相同。

### 6. 从输入端与输出端统一理解

对于

```math
A:\mathbb{R}^n\longrightarrow\mathbb{R}^m,
```

$`N(A)`$ 与 $`\mathrm{Row}(A)`$ 都位于输入空间 $`\mathbb{R}^n`$；左乘 $`P`$ 只改变输出端，所以这两个空间不变。

$`\mathrm{Col}(A)`$ 与 $`N(A^T)`$ 都位于输出空间 $`\mathbb{R}^m`$；左乘 $`P`$ 改变输出坐标，所以这两个空间随 $`P`$ 一起变换。

```text
输入空间 R^n                         输出空间 R^m
N(A)       不变                      Col(A)       变为 P Col(A)
Row(A)     不变                      N(A^T)       变为 P N(A^T)
```

这也解释了为什么可以用 RREF 求原矩阵的秩、零空间和行空间，却不能直接把 R 的主元列当作原矩阵 $`A`$ 的列空间基；行变换已经改变了列空间的位置，求列空间基时必须回到原矩阵取对应列。

## 第 4 题：相同四个基本子空间是否推出 $`A=cB`$

### 1. 判断结论

设 $`A,B\in\mathbb{R}^{m\times n}`$，并且它们的四个基本子空间分别相同：

```math
\mathrm{Col}(A)=\mathrm{Col}(B),
\qquad
N(A)=N(B),
```

```math
\mathrm{Col}(A^T)=\mathrm{Col}(B^T),
\qquad
N(A^T)=N(B^T).
```

判断是否必然存在一个常数 $`c`$，使得

```math
A=cB.
```

结论是：

> **不一定。相同的四个基本子空间不能推出 $`A=cB`$。**

四个子空间记录的是方向、可达性和被压为零的方向；$`A=cB`$ 则要求矩阵的每一个元素都由同一个比例因子统一改变，条件强得多。

### 2. 最简单的可逆矩阵反例

取

```math
A=I_2=\begin{bmatrix}1&0\\0&1\end{bmatrix},
\qquad
B=\begin{bmatrix}1&1\\0&1\end{bmatrix}.
```

两个矩阵都可逆，因此它们的秩都是 $`2`$。对任意可逆 $`2\times2`$ 矩阵，都有

```math
\mathrm{Col}(A)=\mathbb{R}^2,
\qquad
N(A)=\{\mathbf{0}\},
```

```math
\mathrm{Col}(A^T)=\mathbb{R}^2,
\qquad
N(A^T)=\{\mathbf{0}\}.
```

所以 $`A`$ 与 $`B`$ 的四个基本子空间全部相同。

但是它们不可能满足 $`A=cB`$。如果比较左上角元素，会得到 $`c=1`$；此时 $`cB`$ 的右上角元素为 $`1`$，而 $`A`$ 的右上角元素为 $`0`$，矛盾。

因此

> **结论：** 四个基本子空间相同，不能推出 $`A=cB`$。

### 3. 为什么任意两个可逆方阵都能构成这类反例

对任意可逆 $`n\times n`$ 矩阵，$`r=n`$，所以

```math
\mathrm{Col}(A)=\mathbb{R}^n,
\qquad
N(A)=\{\mathbf{0}\}.
```

转置仍然可逆，因此

```math
\mathrm{Col}(A^T)=\mathbb{R}^n,
\qquad
N(A^T)=\{\mathbf{0}\}.
```

也就是说，所有可逆 $`n\times n`$ 矩阵的四个基本子空间都相同，但可逆矩阵有很多，它们显然不可能全都只是彼此的常数倍。

### 4. 非满秩反例

即使不使用可逆矩阵，也可以取

```math
A=\begin{bmatrix}1&0&0\\0&2&0\\0&0&0\end{bmatrix},
\qquad
B=\begin{bmatrix}3&0&0\\0&4&0\\0&0&0\end{bmatrix}.
```

它们的四个基本子空间分别都是

```math
\mathrm{Col}(A)=\mathrm{Col}(B)
=\mathrm{span}\{\mathbf{e}_1,\mathbf{e}_2\},
```

```math
N(A)=N(B)=\mathrm{span}\{\mathbf{e}_3\},
```

```math
\mathrm{Col}(A^T)=\mathrm{Col}(B^T)
=\mathrm{span}\{\mathbf{e}_1,\mathbf{e}_2\},
```

```math
N(A^T)=N(B^T)=\mathrm{span}\{\mathbf{e}_3\}.
```

如果 $`A=cB`$，第一条对角元素要求

```math
1=3c,
```

而第二条对角元素又要求

```math
2=4c.
```

前者给出 $`c=1/3`$，后者给出 $`c=1/2`$，相互矛盾。因此它们也不是常数倍关系。

### 5. 能推出的正确方向

如果 $`A=cB`$ 且 $`c\neq0`$，那么数乘不会改变列空间、零空间、行空间和左零空间：

```math
\mathrm{Col}(cB)=\mathrm{Col}(B),
\qquad
N(cB)=N(B),
```

```math
\mathrm{Col}((cB)^T)=\mathrm{Col}(B^T),
\qquad
N((cB)^T)=N(B^T).
```

所以，若 $`A=cB`$ 且 $`c\neq0`$，那么四个基本子空间相同。

但反方向不成立。

> **四个基本子空间是矩阵的“方向骨架”，不是矩阵的全部数值信息。它们不能告诉我们每个方向具体被放大多少、被剪切多少。**

## 第 5 题：判断一个方阵是否可逆

### 1. 可逆的前提

通常意义下的逆矩阵只对方阵定义。因此首先要确认

```math
A\in\mathbb{R}^{n\times n}.
```

矩形矩阵仍然可以讨论秩、零空间、列独立和方程的解，但一般不谈一个同时满足左右逆的普通 $`A^{-1}`$。

### 2. 可逆矩阵的等价条件

对于 $`n\times n`$ 方阵 $`A`$，下面条件完全等价：

| 观察角度 | 可逆时的条件 |
| --- | --- |
| 秩 | $`\mathrm{rank}(A)=n`$ |
| 主元 | 每一行、每一列都有主元 |
| RREF | $`\mathrm{rref}(A)=I_n`$ |
| 零空间 | $`N(A)=\{\mathbf{0}\}`$ |
| 列 | 所有列线性无关 |
| 列空间 | $`\mathrm{Col}(A)=\mathbb{R}^n`$ |
| 方程 | 对每个 $`\mathbf{b}\in\mathbb{R}^n`$，$`A\mathbf{x}=\mathbf{b}`$ 都有唯一解 |
| 行列式 | $`\det(A)\neq0`$ |
| 逆矩阵 | 存在 $`A^{-1}`$，且 $`A^{-1}A=AA^{-1}=I_n`$ |

因此，只要证明其中任意一条，就可以推出 $`A`$ 可逆。

### 3. 用主元判断

例如

```math
A=\begin{bmatrix}1&2\\3&4\end{bmatrix}
\longrightarrow
\begin{bmatrix}1&2\\0&-2\end{bmatrix}.
```

消元后有两个主元，所以

```math
r=2=n.
```

继续消元得到

```math
\mathrm{rref}(A)=I_2,
```

因此 $`A`$ 可逆。

反例

```math
B=\begin{bmatrix}1&2\\2&4\end{bmatrix}
\longrightarrow
\begin{bmatrix}1&2\\0&0\end{bmatrix}
```

只有一个主元：

```math
r=1<2.
```

所以 $`B`$ 不可逆。原因也可以从列关系看出：

```math
\begin{bmatrix}2\\4\end{bmatrix}
=2\begin{bmatrix}1\\2\end{bmatrix}.
```

列相关意味着存在非零向量

```math
\mathbf{x}=\begin{bmatrix}-2\\1\end{bmatrix}
```

使得

```math
B\mathbf{x}=\mathbf{0}.
```

于是 $`N(B)\neq\{\mathbf{0}\}`$。矩阵把一个非零输入压成了零，信息已经丢失，因此不可能再用逆矩阵唯一恢复输入。

### 4. 一对一与到达全部输出

把方阵看成线性变换

```math
A:\mathbb{R}^n\longrightarrow\mathbb{R}^n.
```

如果

```math
N(A)=\{\mathbf{0}\},
```

则 $`A`$ 是 one-to-one。因为若 $`A\mathbf{x}_1=A\mathbf{x}_2`$，就有

```math
A(\mathbf{x}_1-\mathbf{x}_2)=\mathbf{0},
```

从而 $`\mathbf{x}_1-\mathbf{x}_2=\mathbf{0}`$，即 $`\mathbf{x}_1=\mathbf{x}_2`$。

方阵中 one-to-one 会自动推出 onto，也就是

```math
\mathrm{Col}(A)=\mathbb{R}^n.
```

因此，对方阵来说，列独立既保证解的唯一性，也保证每个右端向量都能到达：

```math
A\text{ 是 }n\times n,
\qquad
r=n
```

推出

```math
N(A)=\{\mathbf{0}\}
\quad\Longrightarrow\quad
\text{解至多一个},
```

```math
\mathrm{Col}(A)=\mathbb{R}^n
\quad\Longrightarrow\quad
\text{每个右端向量都有解}.
```

合起来就是：

```math
\boxed{
A\text{ 为方阵且 }r=n
\quad\Longrightarrow\quad
A\mathbf{x}=\mathbf{b}
\text{ 对每个 }\mathbf{b}\text{ 都有唯一解}
}.
```

### 5. 最实用的判断链

```text
A 是方阵
    -> 检查每一列是否都有 pivot
    -> r=n
    -> N(A)={0}
    -> 列线性无关
    -> Col(A)=R^n
    -> Ax=b 对任意 b 都有唯一解
    -> A 可逆
```

## 第 6 题：正交、点积与四个基本子空间

### 1. $`\mathbf{x}^T\mathbf{y}`$ 就是点积

若 $`\mathbf{x},\mathbf{y}\in\mathbb{R}^n`$ 都写成列向量，直接相乘 $`\mathbf{x}\mathbf{y}`$ 不符合矩阵乘法的尺寸要求。先转置：

```math
\mathbf{x}^T\mathbf{y}
=
\begin{bmatrix}x_1&\cdots&x_n\end{bmatrix}
\begin{bmatrix}y_1\\\vdots\\y_n\end{bmatrix}
=x_1y_1+\cdots+x_ny_n.
```

这正是熟悉的点积：

```math
\boxed{\mathbf{x}\cdot\mathbf{y}=\mathbf{x}^T\mathbf{y}}.
```

因此以下三种说法完全等价：

```math
\mathbf{x}\perp\mathbf{y}
\quad\Longleftrightarrow\quad
\mathbf{x}\cdot\mathbf{y}=0
\quad\Longleftrightarrow\quad
\mathbf{x}^T\mathbf{y}=0.
```

例如

```math
\mathbf{x}=\begin{bmatrix}1\\2\end{bmatrix},
\qquad
\mathbf{y}=\begin{bmatrix}2\\-1\end{bmatrix}
```

满足

```math
\mathbf{x}^T\mathbf{y}=1\cdot2+2\cdot(-1)=0,
```

所以两个向量正交。

### 2. 为什么零空间与行空间正交

把 $`A`$ 的各行记作 $`\mathbf{r}_1^T,\ldots,\mathbf{r}_m^T`$。若

```math
\mathbf{x}\in N(A),
```

那么

```math
A\mathbf{x}=\mathbf{0}
```

按行展开就是

```math
\mathbf{r}_1^T\mathbf{x}=0,
\qquad
\ldots,
\qquad
\mathbf{r}_m^T\mathbf{x}=0.
```

所以 $`\mathbf{x}`$ 与 $`A`$ 的每一行都正交。

而行空间是这些行的所有线性组合：

```math
\mathrm{Row}(A)=\mathrm{Col}(A^T)
=\mathrm{span}\{\mathbf{r}_1,\ldots,\mathbf{r}_m\}.
```

任取行空间中的向量

```math
\mathbf{y}=c_1\mathbf{r}_1+\cdots+c_m\mathbf{r}_m,
```

就有

```math
\mathbf{y}^T\mathbf{x}
=c_1\mathbf{r}_1^T\mathbf{x}+\cdots+c_m\mathbf{r}_m^T\mathbf{x}
=0.
```

因此

```math
\boxed{N(A)\perp\mathrm{Col}(A^T)}.
```

关键逻辑是：

```text
x 属于 N(A)
    -> x 与 A 的每一行正交
    -> x 与这些行的任意线性组合正交
    -> x 与整个 Row(A)=Col(A^T) 正交
```

同理可以得到另一对：

```math
\boxed{\mathrm{Col}(A)\perp N(A^T)}.
```

### 3. 正交子空间的交集

两个子空间 $`S,T`$ 正交，意思是任意 $`\mathbf{s}\in S`$ 和任意 $`\mathbf{t}\in T`$ 都满足

```math
\mathbf{s}^T\mathbf{t}=0.
```

如果某个向量同时属于两者，那么它必须和自己正交：

```math
\mathbf{v}\in S\cap T
\quad\Longrightarrow\quad
\mathbf{v}^T\mathbf{v}=0.
```

但

```math
\mathbf{v}^T\mathbf{v}
=v_1^2+\cdots+v_n^2
```

只有在 $`\mathbf{v}=\mathbf{0}`$ 时才为零。因此

```math
S\perp T
\quad\Longrightarrow\quad
S\cap T=\{\mathbf{0}\}.
```

零向量与任何向量都正交，因为 $`\mathbf{0}^T\mathbf{x}=0`$；所以“正交子空间只能在零向量处相交”并不矛盾。

### 4. 正交补与维数

若在 $`\mathbb{R}^3`$ 中，行空间是一条直线，那么它的正交补必须是一个二维平面；零空间不可能也是一条直线。因为

```math
\dim\bigl(\mathrm{Row}(A)\bigr)=r,
\qquad
\dim N(A)=n-r,
```

所以

```math
\boxed{
\dim\bigl(\mathrm{Row}(A)\bigr)+\dim N(A)=n
}.
```

当 $`n=3`$ 时，可能的维数分配是 $`1+2`$ 或 $`2+1`$，而不是 $`1+1`$。

这也可以写成正交补关系：

```math
N(A)=\mathrm{Row}(A)^\perp,
\qquad
\mathrm{Row}(A)=N(A)^\perp.
```

特殊情况下：

```math
(\mathbb{R}^n)^\perp=\{\mathbf{0}\},
\qquad
\{\mathbf{0}\}^\perp=\mathbb{R}^n.
```

因此一条非零直线不可能与整个 $`\mathbb{R}^2`$ 正交，因为整个空间包含这条直线自身的非零方向；但任何子空间都与零子空间正交。

## 第 7 题：从无解方程到最小二乘正规方程

### 1. 原方程无解时怎么办

因为

```math
A\mathbf{x}\in\mathrm{Col}(A),
```

如果

```math
\mathbf{b}\notin\mathrm{Col}(A),
```

那么

```math
A\mathbf{x}=\mathbf{b}
```

无解。最小二乘法不再要求完全达到 $`\mathbf{b}`$，而是在列空间中寻找距离 $`\mathbf{b}`$ 最近的点：

```math
\mathbf{p}=A\hat{\mathbf{x}}\in\mathrm{Col}(A).
```

令误差为

```math
\mathbf{e}=\mathbf{b}-A\hat{\mathbf{x}}.
```

### 2. 最近点为什么产生正交条件

如果误差还包含沿列空间的方向，就可以在列空间里继续移动 $`A\hat{\mathbf{x}}`$，使它离 $`\mathbf{b}`$ 更近。因此最近点必须满足

```math
\mathbf{e}\perp\mathrm{Col}(A).
```

而

```math
\mathrm{Col}(A)^\perp=N(A^T),
```

所以

```math
\mathbf{e}\in N(A^T)
\quad\Longleftrightarrow\quad
A^T\mathbf{e}=\mathbf{0}.
```

代入误差的定义：

```math
A^T(\mathbf{b}-A\hat{\mathbf{x}})=\mathbf{0}.
```

整理得到

```math
\boxed{A^TA\hat{\mathbf{x}}=A^T\mathbf{b}}.
```

这就是**正规方程**（normal equations）。它不是把原本无解的 $`A\mathbf{x}=\mathbf{b}`$ 变成了完全相等，而是把问题改成了列空间中的正交投影。

### 3. 一个最小例子

取

```math
A=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
\mathbf{b}=\begin{bmatrix}1\\0\end{bmatrix}.
```

此时

```math
\mathrm{Col}(A)
=\mathrm{span}\left\{\begin{bmatrix}1\\1\end{bmatrix}\right\},
```

而 $`\mathbf{b}`$ 不在这条直线上，所以原方程无解。

正规方程中的两项为

```math
A^TA=\begin{bmatrix}1&1\end{bmatrix}\begin{bmatrix}1\\1\end{bmatrix}=2,
\qquad
A^T\mathbf{b}=\begin{bmatrix}1&1\end{bmatrix}\begin{bmatrix}1\\0\end{bmatrix}=1.
```

于是

```math
2\hat{x}=1,
\qquad
\hat{x}=\frac12.
```

投影点为

```math
A\hat{x}
=\begin{bmatrix}1/2\\1/2\end{bmatrix},
```

误差为

```math
\mathbf{e}
=\begin{bmatrix}1\\0\end{bmatrix}
-\begin{bmatrix}1/2\\1/2\end{bmatrix}
=\begin{bmatrix}1/2\\-1/2\end{bmatrix}.
```

检查正交性：

```math
\begin{bmatrix}1&1\end{bmatrix}
\begin{bmatrix}1/2\\-1/2\end{bmatrix}=0.
```

所以误差确实垂直于列空间。

### 4. 正规方程的逻辑链

```text
Ax=b 无解
    -> b 不在 Col(A) 中
    -> 在 Col(A) 中找最近点 A x_hat
    -> 误差 e=b-A x_hat 垂直于 Col(A)
    -> e 属于 N(A^T)
    -> A^T e=0
    -> A^T A x_hat=A^T b
```

如果 $`A`$ 的列线性无关，正规方程有唯一解；这需要下一题的等价关系。

## 第 8 题：为什么 $`A^TA`$ 可逆等价于列独立

### 1. 先确认尺寸

若

```math
A\in\mathbb{R}^{m\times n},
```

那么

```math
A^TA\in\mathbb{R}^{n\times n}.
```

所以 $`A^TA`$ 是方阵，可以使用“方阵可逆当且仅当零空间只有零向量”的判断。

### 2. 证明两个零空间相同

先假设

```math
A^TA\mathbf{x}=\mathbf{0}.
```

左乘 $`\mathbf{x}^T`$：

```math
\mathbf{x}^TA^TA\mathbf{x}=0.
```

利用

```math
\mathbf{x}^TA^T=(A\mathbf{x})^T,
```

可得

```math
(A\mathbf{x})^T(A\mathbf{x})=0.
```

而向量与自身的点积是长度平方：

```math
\|A\mathbf{x}\|^2=0.
```

实向量长度平方为零只能说明

```math
A\mathbf{x}=\mathbf{0}.
```

反过来，如果 $`A\mathbf{x}=\mathbf{0}`$，显然

```math
A^TA\mathbf{x}=A^T\mathbf{0}=\mathbf{0}.
```

因此

```math
\boxed{A^TA\mathbf{x}=\mathbf{0}
\quad\Longleftrightarrow\quad
A\mathbf{x}=\mathbf{0}}.
```

也就是

```math
\boxed{N(A^TA)=N(A)}.
```

### 3. 接回列线性无关

前面已经知道

```math
\text{A 的列线性无关}
\quad\Longleftrightarrow\quad
N(A)=\{\mathbf{0}\}.
```

结合 $`N(A^TA)=N(A)`$，得到

```math
\begin{aligned}
A\text{ 的列线性无关}
&\Longleftrightarrow N(A)=\{\mathbf{0}\}\\
&\Longleftrightarrow N(A^TA)=\{\mathbf{0}\}\\
&\Longleftrightarrow A^TA\text{ 可逆}.
\end{aligned}
```

最终结论是

```math
\boxed{
A^TA\text{ 可逆}
\quad\Longleftrightarrow\quad
A\text{ 的列向量线性无关}
}.
```

### 4. 为什么这对最小二乘重要

正规方程是

```math
A^TA\hat{\mathbf{x}}=A^T\mathbf{b}.
```

如果 $`A`$ 的列独立，那么 $`A^TA`$ 可逆，于是可以唯一求出

```math
\boxed{
\hat{\mathbf{x}}
=(A^TA)^{-1}A^T\mathbf{b}
}.
```

如果 $`A`$ 的列相关，那么 $`N(A)\neq\{\mathbf{0}\}`$，从而 $`N(A^TA)\neq\{\mathbf{0}\}`$，$`A^TA`$ 不可逆，不能直接写出上面的逆矩阵公式。

这条结论把前面的列独立、零空间、可逆矩阵和最小二乘解连接起来：

```text
A 的列独立
    -> N(A)={0}
    -> N(A^T A)={0}
    -> A^T A 可逆
    -> 正规方程有唯一解
    -> 最小二乘解可以写成 (A^T A)^(-1) A^T b
```

## Review 1 当前总结

目前的习题内容从八个方向使用了前面学过的零空间、列空间、秩、正交与可逆性：

1. 第 1 题从“特解 + 零空间”反推矩阵，并由列空间判断哪些右端向量可解。
2. 第 2 题利用可逆左乘不改变零空间，把复杂矩阵的齐次方程转化为更简单矩阵的齐次方程。
3. 第 3 题区分行变换后保持不变的输入端空间，以及随输出坐标变化的输出端空间。
4. 第 4 题说明四个基本子空间只描述结构，不足以确定矩阵的具体数值，更不能推出两个矩阵是常数倍。
5. 第 5 题把方阵可逆与满秩、主元、零空间、列独立、`Ax=b` 的存在唯一性以及行列式统一起来。
6. 第 6 题从 $`x^T y`$ 的点积形式出发，证明零空间与行空间、列空间与左零空间之间的正交关系。
7. 第 7 题从 $`Ax=b`$ 无解出发，用误差向量垂直于列空间的条件推出最小二乘正规方程。
8. 第 8 题证明 $`N(A^T A)=N(A)`$，从而说明 $`A^T A`$ 可逆当且仅当 $`A`$ 的列向量线性无关。

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
N(A)\perp\mathrm{Col}(A^T),
\qquad
\mathrm{Col}(A)\perp N(A^T),
```

```math
A^T A\widehat{\mathbf{x}}=A^T\mathbf{b},
\qquad
N(A^T A)=N(A).
```

对于方阵，判断可逆可以沿着下面的等价链理解：

```text
A 是 n×n 方阵
    -> r=n
    -> 每行、每列都有 pivot
    -> rref(A)=I_n
    -> N(A)={0}，列向量线性无关
    -> Col(A)=R^n
    -> 对每个 b 都有唯一解 Ax=b
    -> A 可逆
```

> 当 $`C`$ 可逆时，$`N(CD)=N(D)`$。

以及行变换 $`B=PA`$ 下的关系

```math
N(B)=N(A),
\qquad
\mathrm{Row}(B)=\mathrm{Row}(A),
```

```math
\mathrm{Col}(B)=P\,\mathrm{Col}(A),
\qquad
N(B^T)=P\,N(A^T).
```

> **Review 1 的主线不是孤立地计算矩阵，而是从解的结构识别零空间，从零空间识别自由度与列关系，再用列空间、秩、行变换、正交和可逆性判断方程的存在性、唯一性、近似解以及四个基本子空间的关系。**
