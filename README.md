# Matrix Toolbox in C

A console-based **Matrix Toolbox** written in C.  
It allows you to create, manipulate, and analyze matrices with essential linear algebra operations.  
The program uses a simple menu system and supports both manual input and random generation.

---

## Features

- Create matrices:
  - Manual input
  - Random initialization
  - Load from file (`.txt`)
- Save matrix to file
- Basic operations:
  - Addition: \( C = A + B \)
  - Subtraction: \( C = A - B \)
  - Multiplication: \( C = A \times B \)
  - Transpose: \( A^T \)
- Determinant:
  \[
  \det(A) = \prod_{i=1}^{n} a_{ii}
  \]
  (computed via Gaussian elimination)
- Inverse matrix (Gauss-Jordan elimination)
- Memory-safe with dynamic allocation

---

## 🛠 Compilation & Run

```bash
gcc -std=c11 -O2 -Wall -lm -o matrix matrix.c
./matrix
````

---

## File Format

Text file representation:

```
rows cols
a11 a12 ... a1n
a21 a22 ... a2n
...
am1 am2 ... amn
```

Example:

```
2 2
1 2
3 4
```

---

## Example Workflow

1. Start program:

   ```
   ./matrix
   ```

2. Create a random $3 \times 3$ matrix.

   Example matrix $A$:

   $$
   A =
   \begin{bmatrix}
   1 & 2 & 3 \\
   4 & 5 & 6 \\
   7 & 8 & 9
   \end{bmatrix}
   $$

3. Compute determinant:

   $$
   \det(A) = 0
   $$

   (matrix is singular, no inverse exists)

4. Try with another matrix:

   $$
   B =
   \begin{bmatrix}
   2 & 1 \\
   5 & 3
   \end{bmatrix}
   $$

   Determinant:

   $$
   \det(B) = (2 \cdot 3) - (5 \cdot 1) = 1
   $$

   Inverse:

   $$
   B^{-1} =
   \begin{bmatrix}
   3 & -1 \\
   -5 & 2
   \end{bmatrix}
   $$

---

## Concepts

* **Transpose**: Flip over diagonal, $(A^T)_{ij} = A_{ji}$.
* **Determinant**: Scalar value showing invertibility. If $\det(A) = 0$, matrix is singular.
* **Inverse**: Matrix $A^{-1}$ such that:

  $$
  A \cdot A^{-1} = I
  $$

---

## Notes

* Only square matrices have determinants and inverses.
* Random matrices are generated within user-defined range.
* Gaussian elimination is $O(n^3)$; efficient for small/medium matrices.

---

## Future Improvements

* Eigenvalues & eigenvectors
* LU decomposition
* Binary file format for faster I/O
