/* matrix.c
   Программа для работы с матрицами — однофайловая.
   Поддерживает: ввод, случайная генерация, вывод, сложение, вычитание,
   умножение, транспонирование, детерминант (через Gaussian elimination),
   обратная матрица (Gauss-Jordan), сохранение/загрузка.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define EPS 1e-12

typedef struct {
    size_t rows;
    size_t cols;
    double *data; // contiguous storage: data[i*cols + j]
} Matrix;

/* ====== Вспомогательные функции для работы с матрицами ====== */

Matrix *matrix_create(size_t rows, size_t cols) {
    Matrix *m = malloc(sizeof(Matrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->data = calloc(rows * cols, sizeof(double));
    if (!m->data) { free(m); return NULL; }
    return m;
}

void matrix_free(Matrix *m) {
    if (!m) return;
    free(m->data);
    free(m);
}

double matrix_get(const Matrix *m, size_t i, size_t j) {
    return m->data[i * m->cols + j];
}

void matrix_set(Matrix *m, size_t i, size_t j, double v) {
    m->data[i * m->cols + j] = v;
}

void matrix_print(const Matrix *m) {
    if (!m) { printf("(null)\n"); return; }
    printf("Matrix %zux%zu:\n", m->rows, m->cols);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            printf("%10.4g ", matrix_get(m, i, j));
        }
        printf("\n");
    }
}

/* Ввод матрицы вручную */
void matrix_input(Matrix *m) {
    printf("Ввод матрицы %zux%zu (по элементам):\n", m->rows, m->cols);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            double v;
            printf("A[%zu][%zu] = ", i, j);
            while (scanf("%lf", &v) != 1) {
                while (getchar() != '\n'); // чистим ввод
                printf("Неверный ввод. Попробуйте снова: ");
            }
            matrix_set(m, i, j, v);
        }
    }
}

/* Заполнение случайными числами в диапазоне [minv, maxv] */
void matrix_random(Matrix *m, double minv, double maxv) {
    for (size_t i = 0; i < m->rows; ++i)
        for (size_t j = 0; j < m->cols; ++j)
            matrix_set(m, i, j, minv + (maxv - minv) * (rand() / (double)RAND_MAX));
}

/* Копирование */
Matrix *matrix_clone(const Matrix *a) {
    Matrix *b = matrix_create(a->rows, a->cols);
    if (!b) return NULL;
    memcpy(b->data, a->data, sizeof(double) * a->rows * a->cols);
    return b;
}

/* Сложение/вычитание */
Matrix *matrix_add_sub(const Matrix *a, const Matrix *b, int subtract) {
    if (!a || !b) return NULL;
    if (a->rows != b->rows || a->cols != b->cols) return NULL;
    Matrix *c = matrix_create(a->rows, a->cols);
    if (!c) return NULL;
    for (size_t i = 0; i < a->rows * a->cols; ++i)
        c->data[i] = a->data[i] + (subtract ? -b->data[i] : b->data[i]);
    return c;
}

/* Умножение */
Matrix *matrix_multiply(const Matrix *a, const Matrix *b) {
    if (!a || !b) return NULL;
    if (a->cols != b->rows) return NULL;
    Matrix *c = matrix_create(a->rows, b->cols);
    if (!c) return NULL;
    for (size_t i = 0; i < a->rows; ++i) {
        for (size_t k = 0; k < a->cols; ++k) {
            double aik = matrix_get(a, i, k);
            for (size_t j = 0; j < b->cols; ++j) {
                c->data[i * c->cols + j] += aik * matrix_get(b, k, j);
            }
        }
    }
    return c;
}

/* Транспонирование */
Matrix *matrix_transpose(const Matrix *a) {
    Matrix *t = matrix_create(a->cols, a->rows);
    if (!t) return NULL;
    for (size_t i = 0; i < a->rows; ++i)
        for (size_t j = 0; j < a->cols; ++j)
            matrix_set(t, j, i, matrix_get(a, i, j));
    return t;
}

/* Сохранение/загрузка в простой текстовый формат:
   Первая строка: rows cols
   Далее rows строк по cols чисел.
*/
int matrix_save_txt(const Matrix *m, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) return 0;
    fprintf(f, "%zu %zu\n", m->rows, m->cols);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            fprintf(f, "%.12g ", matrix_get(m, i, j));
        }
        fprintf(f, "\n");
    }
    fclose(f);
    return 1;
}

Matrix *matrix_load_txt(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) return NULL;
    size_t rows, cols;
    if (fscanf(f, "%zu %zu", &rows, &cols) != 2) { fclose(f); return NULL; }
    Matrix *m = matrix_create(rows, cols);
    if (!m) { fclose(f); return NULL; }
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            if (fscanf(f, "%lf", &m->data[i * cols + j]) != 1) {
                matrix_free(m);
                fclose(f);
                return NULL;
            }
    fclose(f);
    return m;
}

/* ====== Линейная алгебра: детерминант и обратная матрица ====== */

/* Вычисление детерминанта квадратной матрицы методом приведения к верхней треугольной форме.
   Возвращает 0 если не квадратная. Работает с копией матрицы (не изменяет входную).
*/
double matrix_determinant(const Matrix *a) {
    if (!a) return 0.0;
    if (a->rows != a->cols) {
        fprintf(stderr, "Determinant: matrix is not square\n");
        return 0.0;
    }
    size_t n = a->rows;
    // Копируем в рабочую матрицу
    double *mat = malloc(n * n * sizeof(double));
    if (!mat) return 0.0;
    memcpy(mat, a->data, n * n * sizeof(double));

    double det = 1.0;
    int sign = 1;
    for (size_t i = 0; i < n; ++i) {
        // Поиск опорного элемента (pivot)
        size_t piv = i;
        for (size_t r = i; r < n; ++r) {
            if (fabs(mat[r*n + i]) > fabs(mat[piv*n + i])) piv = r;
        }
        if (fabs(mat[piv*n + i]) < EPS) {
            det = 0.0;
            break;
        }
        if (piv != i) {
            // swap rows i and piv
            for (size_t c = 0; c < n; ++c) {
                double tmp = mat[i*n + c];
                mat[i*n + c] = mat[piv*n + c];
                mat[piv*n + c] = tmp;
            }
            sign = -sign;
        }
        double pivot = mat[i*n + i];
        det *= pivot;
        // нормируем и вычитаем
        for (size_t r = i + 1; r < n; ++r) {
            double factor = mat[r*n + i] / pivot;
            for (size_t c = i; c < n; ++c) {
                mat[r*n + c] -= factor * mat[i*n + c];
            }
        }
    }
    det *= sign;
    free(mat);
    return det;
}

/* Обратная матрица методом Гаусса-Жордана.
   Возвращает NULL, если матрица не квадратная или необратима.
*/
Matrix *matrix_inverse(const Matrix *a) {
    if (!a) return NULL;
    if (a->rows != a->cols) {
        fprintf(stderr, "Inverse: matrix is not square\n");
        return NULL;
    }
    size_t n = a->rows;
    // Создаём расширенную матрицу nx(2n)
    double *E = malloc(n * 2 * n * sizeof(double));
    if (!E) return NULL;
    // Заполняем [A | I]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            E[i*(2*n) + j] = matrix_get(a, i, j);
            E[i*(2*n) + (n + j)] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Прямой ход
    for (size_t i = 0; i < n; ++i) {
        // pivot selection
        size_t piv = i;
        for (size_t r = i; r < n; ++r)
            if (fabs(E[r*(2*n) + i]) > fabs(E[piv*(2*n) + i])) piv = r;
        if (fabs(E[piv*(2*n) + i]) < EPS) {
            free(E);
            return NULL; // сингулярная матрица
        }
        // swap rows if needed
        if (piv != i) {
            for (size_t c = 0; c < 2*n; ++c) {
                double tmp = E[i*(2*n) + c];
                E[i*(2*n) + c] = E[piv*(2*n) + c];
                E[piv*(2*n) + c] = tmp;
            }
        }
        // normalize row i
        double div = E[i*(2*n) + i];
        for (size_t c = 0; c < 2*n; ++c) E[i*(2*n) + c] /= div;
        // eliminate other rows
        for (size_t r = 0; r < n; ++r) {
            if (r == i) continue;
            double factor = E[r*(2*n) + i];
            if (fabs(factor) < EPS) continue;
            for (size_t c = 0; c < 2*n; ++c)
                E[r*(2*n) + c] -= factor * E[i*(2*n) + c];
        }
    }

    // Результат в правой половине
    Matrix *inv = matrix_create(n, n);
    if (!inv) { free(E); return NULL; }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            inv->data[i*n + j] = E[i*(2*n) + (n + j)];

    free(E);
    return inv;
}

/* ====== Меню и взаимодействие с пользователем ====== */

void flush_stdin(void) {
    int c;
    while ((c = getchar()) != '\n' && c != EOF) {}
}

void print_menu(void) {
    puts("\n=== Matrix Toolbox ===");
    puts("1) Создать новую матрицу вручную");
    puts("2) Создать новую матрицу случайно");
    puts("3) Загрузить матрицу из файла");
    puts("4) Показать текущую матрицу");
    puts("5) Сохранить текущую матрицу в файл");
    puts("6) Сложить с другой матрицей");
    puts("7) Вычесть другую матрицу");
    puts("8) Умножить на другую матрицу");
    puts("9) Транспонировать текущую матрицу");
    puts("10) Детерминант (если квадратная)");
    puts("11) Обратная матрица (если квадратная и невырождена)");
    puts("12) Освободить текущую матрицу");
    puts("0) Выход");
    printf("Выберите действие: ");
}

Matrix *ask_create_manual(void) {
    size_t r, c;
    printf("Введите число строк: ");
    while (scanf("%zu", &r) != 1) { flush_stdin(); printf("Неверно. Введите число строк: "); }
    printf("Введите число столбцов: ");
    while (scanf("%zu", &c) != 1) { flush_stdin(); printf("Неверно. Введите число столбцов: "); }
    Matrix *m = matrix_create(r, c);
    if (!m) { fprintf(stderr, "Не удалось выделить память\n"); return NULL; }
    matrix_input(m);
    return m;
}

Matrix *ask_create_random(void) {
    size_t r, c;
    double minv, maxv;
    printf("Введите число строк: ");
    while (scanf("%zu", &r) != 1) { flush_stdin(); printf("Неверно. Введите число строк: "); }
    printf("Введите число столбцов: ");
    while (scanf("%zu", &c) != 1) { flush_stdin(); printf("Неверно. Введите число столбцов: "); }
    printf("Минимум для случайных: ");
    while (scanf("%lf", &minv) != 1) { flush_stdin(); printf("Неверно. Введите число: "); }
    printf("Максимум для случайных: ");
    while (scanf("%lf", &maxv) != 1) { flush_stdin(); printf("Неверно. Введите число: "); }
    if (maxv < minv) { double t = minv; minv = maxv; maxv = t; }
    Matrix *m = matrix_create(r, c);
    if (!m) { fprintf(stderr, "Не удалось выделить память\n"); return NULL; }
    matrix_random(m, minv, maxv);
    return m;
}

Matrix *ask_load_file(void) {
    char fname[512];
    printf("Имя файла для загрузки: ");
    scanf("%511s", fname);
    Matrix *m = matrix_load_txt(fname);
    if (!m) fprintf(stderr, "Не удалось загрузить матрицу из '%s'\n", fname);
    return m;
}

int ask_save_file(const Matrix *m) {
    char fname[512];
    printf("Имя файла для сохранения: ");
    scanf("%511s", fname);
    if (matrix_save_txt(m, fname)) {
        printf("Сохранено в '%s'\n", fname);
        return 1;
    } else {
        fprintf(stderr, "Ошибка при сохранении в '%s'\n", fname);
        return 0;
    }
}

Matrix *ask_other_matrix_for_operation(void) {
    puts("Выберите способ задания второй матрицы:");
    puts("1) Ввести вручную");
    puts("2) Сгенерировать случайно");
    puts("3) Загрузить из файла");
    printf("Выбор: ");
    int choice;
    if (scanf("%d", &choice) != 1) { flush_stdin(); return NULL; }
    if (choice == 1) return ask_create_manual();
    if (choice == 2) return ask_create_random();
    if (choice == 3) return ask_load_file();
    return NULL;
}

int main(void) {
    srand((unsigned)time(NULL));
    Matrix *M = NULL;
    int running = 1;
    while (running) {
        print_menu();
        int opt;
        if (scanf("%d", &opt) != 1) { flush_stdin(); continue; }
        switch (opt) {
            case 1:
                if (M) { matrix_free(M); M = NULL; }
                M = ask_create_manual();
                break;
            case 2:
                if (M) { matrix_free(M); M = NULL; }
                M = ask_create_random();
                break;
            case 3:
                if (M) { matrix_free(M); M = NULL; }
                M = ask_load_file();
                break;
            case 4:
                if (!M) printf("Текущая матрица отсутствует.\n");
                else matrix_print(M);
                break;
            case 5:
                if (!M) { printf("Нет матрицы для сохранения.\n"); break; }
                ask_save_file(M);
                break;
            case 6: { // add
                if (!M) { printf("Нет текущей матрицы.\n"); break; }
                Matrix *B = ask_other_matrix_for_operation();
                if (!B) { printf("Операция отменена.\n"); break; }
                Matrix *C = matrix_add_sub(M, B, 0);
                if (!C) printf("Ошибка: несовместимые размеры или память.\n");
                else { printf("Результат (сложение):\n"); matrix_print(C); matrix_free(C); }
                matrix_free(B);
                break;
            }
            case 7: { // sub
                if (!M) { printf("Нет текущей матрицы.\n"); break; }
                Matrix *B = ask_other_matrix_for_operation();
                if (!B) { printf("Операция отменена.\n"); break; }
                Matrix *C = matrix_add_sub(M, B, 1);
                if (!C) printf("Ошибка: несовместимые размеры или память.\n");
                else { printf("Результат (вычитание):\n"); matrix_print(C); matrix_free(C); }
                matrix_free(B);
                break;
            }
            case 8: { // mul
                if (!M) { printf("Нет текущей матрицы.\n"); break; }
                Matrix *B = ask_other_matrix_for_operation();
                if (!B) { printf("Операция отменена.\n"); break; }
                Matrix *C = matrix_multiply(M, B);
                if (!C) printf("Ошибка: несовместимые размеры или память.\n");
                else { printf("Результат (умножение):\n"); matrix_print(C); matrix_free(C); }
                matrix_free(B);
                break;
            }
            case 9: { // transpose
                if (!M) { printf("Нет текущей матрицы.\n"); break; }
                Matrix *T = matrix_transpose(M);
                if (!T) printf("Ошибка: память.\n");
                else {
                    matrix_free(M);
                    M = T;
                    printf("Транспонирование выполнено. Теперь матрица имеет размер %zux%zu\n", M->rows, M->cols);
                }
                break;
            }
            case 10: { // determinant
                if (!M) { printf("Нет текущей матрицы.\n"); break; }
                if (M->rows != M->cols) { printf("Не квадратная матрица.\n"); break; }
                double det = matrix_determinant(M);
                printf("Детерминант = %.12g\n", det);
                break;
            }
            case 11: { // inverse
                if (!M) { printf("Нет текущей матрицы.\n"); break; }
                if (M->rows != M->cols) { printf("Не квадратная матрица.\n"); break; }
                Matrix *inv = matrix_inverse(M);
                if (!inv) printf("Матрица необратима или ошибка.\n");
                else { printf("Обратная матрица:\n"); matrix_print(inv); matrix_free(inv); }
                break;
            }
            case 12:
                if (M) { matrix_free(M); M = NULL; printf("Матрица освобождена.\n"); }
                else printf("Матрица отсутствует.\n");
                break;
            case 0:
                running = 0;
                break;
            default:
                printf("Неизвестный пункт меню.\n");
        }
    }

    if (M) matrix_free(M);
    puts("Выход. Пока!");
    return 0;
}

