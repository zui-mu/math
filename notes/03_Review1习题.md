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

