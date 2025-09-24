# Matrix Toolbox

A C console program for working with matrices.  
Supports: creation (manual/random), file I/O, addition, subtraction, multiplication, transpose, determinant, and inverse.

---

## Mathematical Background

### Matrix Example (3×3)
\[
A =
\begin{bmatrix}
1 & 2 & 3 \\
4 & 5 & 6 \\
7 & 8 & 9
\end{bmatrix}
\]

Determinant:
\[
\det(A) = 0
\]
(Singular → no inverse exists)

---

### Matrix Example (2×2)
\[
B =
\begin{bmatrix}
2 & 1 \\
5 & 3
\end{bmatrix}
\]

Determinant:
\[
\det(B) = (2 \cdot 3) - (5 \cdot 1) = 1
\]

Inverse:
\[
B^{-1} =
\begin{bmatrix}
3 & -1 \\
-5 & 2
\end{bmatrix}
\]

---

## Implemented Operations

- **Addition / Subtraction**  
  \[
  (A \pm B)_{ij} = A_{ij} \pm B_{ij}
  \]

- **Multiplication**  
  \[
  (AB)_{ij} = \sum_{k=1}^{n} A_{ik}B_{kj}
  \]

- **Transpose**  
  \[
  (A^T)_{ij} = A_{ji}
  \]

- **Determinant & Inverse**  
  via Gaussian elimination with partial pivoting.

---

## Complexity

- Addition/Subtraction: \(O(n^2)\)  
- Multiplication: \(O(n^3)\)  
- Determinant/Inverse: \(O(n^3)\)

---

## Usage

Compile:
```bash
gcc -std=c11 -O2 -Wall -lm -o matrix matrix.c
````

Run:

```bash
./matrix
```

Console demo:

```
=== Matrix Toolbox ===
1) Create new matrix
2) Random matrix
...
> 2
Rows: 3
Cols: 3
Min: 0
Max: 10
Random matrix created:
  1.00  2.00  3.00
  4.00  5.00  6.00
  7.00  8.00  9.00

Determinant = 0
```
