#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int n;
int it;
const char *FORM = "%*.5lf";
enum 
{
    FORM_L = 12
};
//печать матрицы
void print_m(double **matrix, char * str) {
    if (str) {
        printf("%s\n", str);
    }
    for (int q = 0; q < n; q++) {
        for (int p = 0; p < n + 1; p ++){
            printf(FORM, FORM_L, matrix[p][q]);
        }    
        printf("\n");
    }
}
//считываение матрицы
void get_m(double **matrix, FILE * f) {
    for (int q = 0; q < n; q++) {
        for (int p = 0; p < n + 1; p ++){
            double a;
            fscanf(f, "%lf", &a);
            matrix[p][q] = a;
        }
    }
}
//копирование матрицы
void copy_m(double **from, double **to) {
    for (int q = 0; q < n; q++) {
        for (int p = 0; p < n + 1; p ++){
            to[p][q] = from[p][q];
        }
    }
}
//выделение памяти под матрицу
double **init_m() {
    double **matrix = (double **)calloc(n + 1, sizeof(matrix));
    for (int t = 0; t < n + 1; t++) {
        matrix[t] = (double *) calloc(n, sizeof(matrix[t]));
    }
    return matrix;
}
//если первый элемент строки 0
int code_zero(double **matrix, double **ed, int t) {
    int j = t + 1;
    while (j < n && (matrix[t][j] == 0)) { //ищем строку где элемент не ноль
        j++;
    }
    if (matrix[t][j] == 0) {
        return 0; // таких нет
    } else {
        for (int l = 0; l < n + 1; l++) { // обмен
            int a = matrix[l][t];
            matrix[l][t] = matrix[l][j];
            matrix[l][j] = a;
            if (ed) {
                a = ed[l][t];
                ed[l][t] = ed[l][j];
                ed[l][j] = a;
            }
            
        }
    }
    return 1;
}
//поиск максимального элемента в строке
void find_max(double **matrix, int *order, int t) {
    double max = fabs(matrix[t][t]);
    int maxn = t;
    for (int p = t + 1; p < n; p++) {
        if (max < fabs(matrix[p][t])) {
            max = fabs(matrix[p][t]);
            maxn = p;
        }
    }
    double * ch = matrix[maxn];
    matrix[maxn] = matrix[t];
    matrix[t] = ch;
    int a = order[maxn];
    order[maxn] = order[t];
    order[t] = a;
}
//инициализация единичной матрицы
void ed_init(double **ed) {
    for(int p = 0; p < n; p++) {
        for (int q = 0; q < n; q++) {
            if (q == p) {
                ed[p][q] = 1;
            } else {
                ed[p][q] = 0;
            }
        }
    }
}
//короткая функция для вычисление определителя
double find_det(double **matrix){
    int det = 1;
    for (int t = 0; t < n; t++) {
        det *= matrix[t][t];
    }
    return det;
}
//прямой метод Гаусса
double str_gauss(double **matrix, double **ed, int mode, int *order) { 
    // mode 0 - нормальный
    // mode 1 - улучшенный
    // mode 2 - обратная матрица
    double det = 1.0;
    for (int t = 0; t < n; t++) {
        if (mode == 1) {
            find_max(matrix, order, t); //ищем максимальный столбец и меняем местами
        }
        if (matrix[t][t] == 0) { // ведущий элемент ноль
            if (code_zero(matrix, ed, t)) {
                det *= -1;
            }   
        }
        double div_el = matrix[t][t]; // делим все на ведущий
        det *= div_el;
        for (int l = t; l < n + 1; l++) {
            matrix[l][t] /= div_el;

        }

        if(mode == 2) { // поиск обратной
            for(int l = 0; l < n; l++){
                ed[l][t] /= div_el;
            }
        }
        // вычитаем строки
        for (int q = t + 1; q < n; q++) {
            double mul = matrix[t][q]; // умножаем на начальный элемент строки
            
            for (int p = t; p < n + 1; p++) {
                matrix[p][q] -= mul * matrix[p][t];
            }
            if (mode == 2) {
                for (int p = 0; p < n + 1; p++) {
                    ed[p][q] -= mul * ed[p][t];
                }
            }
        }
    }
    
    return det * find_det(matrix);
}
//обратный метод Гаусса
void rev_gauss(double **matrix, double **ed, int mode) {
    for (int t = n - 1; t > 0; t--) { // идем по строкам
        for (int q = t - 1; q >= 0; q--) { // вычитаем строки
            double mul = matrix[t][q]; // умножаем на начальный элемент строки
            for (int p = t; p < n + 1; p++) {
                matrix[p][q] -= mul * matrix[p][t];
            }
            if (mode == 2) {
                for (int p = 0; p < n; p++) {
                    ed[p][q] -= mul * ed[p][t];
                }
            }
        }
    }
    
}
//генерирует матрицу по формулам
void gen_m(double **matrix, double x) {
    double qM = 1.001 - 2 * 6 * 0.001;
    for (int q = 0; q < n; q++) {
        for (int p = 0; p < n; p++) {
            if (q != p) {
                matrix[p][q] = pow(qM, p + q + 2) + 0.1 * (p - q);
            } else {
                matrix[p][q] = pow(qM - 1, p + q + 2);
            }
        }
    }
    for (int t = 0; t < n; t++) {
        matrix[n][t] = x * exp(x / (t + 1)) * cos(x / (t + 1));
    }
}
//получает вектор столбец ответов
double *get_res(double **matrix, int *order) {
    double *f = calloc(n, sizeof(f));

    for(int t = 0; t < n; t++){
        if (order) {
            for (int j = 0; j < n; j++){
                if (order[j] == t) {
                    f[j] = matrix[n][t];
                }
            }
        } else {
            f[t] = matrix[n][t];
        }
        
    }
    return f;
}
//печатает массив строку
void print_v(double *f, char *s) {
    if (s) {
        printf("%s\n", s);
    }
    for (int t = 0; t < n; t++) {
        printf(FORM, FORM_L, f[t]);
    }
    printf("\n");
}
//инициализация массива обменов для улучшенного метода Гаусса
int *init_order() {
    int *order = calloc(n, sizeof(n));
    for (int p = 0; p < n; p++){
        order[p] = p;
    }
    return order;
}
//инициализация матрицы для функций оберток
double **init_all(double **m) {
    double **matrix = init_m();
    copy_m(m, matrix);
    print_m(matrix, "matrix:");
    return matrix;
}
//метод Гаусса(обычный)
double **orig_g(double **m) {
    printf("SIMPLE GAUSS METHOD\n");
    double **matrix = init_all(m);

    str_gauss(matrix, NULL, 0, NULL);
    rev_gauss(matrix, NULL, 0);
    
    print_m(matrix, "changed matrix:");
    double *res = get_res(matrix, NULL);
    print_v(res, "answer:");
    printf("\n");
    return matrix;
}
//метод Гаусса(улучшенный)
double **upgr_g(double **m) {
    printf("UPGRADE GAUSS METHOD\n");
    double **matrix = init_all(m);

    int *order = init_order();
    str_gauss(matrix, NULL, 1, order);
    rev_gauss(matrix, NULL, 0);

    print_m(matrix, "changed matrix:");
    double *res = get_res(matrix, order);
    print_v(res, "answer:");
    printf("\n");
    return matrix;   
}
//метод получающий обратную матрицу
double **get_rev(double **m) {
    printf("GETTING REVERSE MATRIX\n");
    double **matrix = init_all(m);
    double **ed;
    ed = init_m();
    ed_init(ed);
    str_gauss(matrix, ed, 2, NULL);
    rev_gauss(matrix, ed, 2);

    print_m(ed, "changed matrix:");
    printf("\n");
    return ed;
}
//освобождение памяти под матрицу
void free_m(double **matrix) {
    for (int t = 0; t < n; t++) {
        free(matrix[t]);
    }
    free(matrix);
}
//вычисление определителя
double get_determinant(double **m) {
    printf("GETTING DETERMINANT\n");
    double **matrix = init_all(m);
    double det = str_gauss(matrix, NULL, 0, NULL);
    printf("det: %lf \n", det);
    free_m(matrix);
    printf("\n");
    return det;
}
//инициализация строки нулями
void string_init_z(double *str){
    for (int t = 0; t < n; t++) {
        str[t] = 0;
    }
}
//вычисляет погрешность МВР
double get_diff(double **matrix, double *res) {
    double sqr = 0.0;
    double r[n];
    for(int t = 0; t < n; t++) {
        r[t] = matrix[n][t];
    }
    for (int q = 0; q < n; q++) {
        for (int p = 0; p < n; p++) {
            r[q] -= matrix[p][q] *  res[p];
        }
        sqr += r[q] * r[q];
    }
    return sqrt(sqr);
}
//вычисляет коэфициент в формуле для МВР
double get_koef(double **matrix, double *pred, double *now,  int q) {

    double koef = matrix[n][q];
    for (int j = 0; j < q; j++) {
        koef -= matrix[j][q] * now[j];
    }
    for (int j = q; j < n; j++) {
        koef -= matrix[j][q] * pred[j];
    }
    return koef;
}
//метод верхней релаксации (основное тело)
double *iter(double **matrix, double om, double eps) {
    double *x = calloc(n, sizeof(x));
    
    double pred[n];
    int i = 0;
    while(i < 1000) {
        for(int q = 0; q < n; q++) {
            double koef;
            for(int t = 0; t < n; t++) {
                pred[t] = x[t];
            }
            koef = get_koef(matrix, pred, x, q);

            x[q] += om * koef / matrix[q][q];
        }
        if(get_diff(matrix, x) < eps){
            printf("w=%*.5lf, iterations =  %*d\n", FORM_L, om, FORM_L, i);
            it = i;
            return x;
        }
        i++;
    }
    return x;
}
//метод верхней релаксации (обертка)
double *it_meth(double **m, double eps) {
    printf("APPLYING ITERATIVE METHOD\n");

    double **matrix; 
    double *x;
    int arr[20];
    double ome[20];
    int j = 0;
    for (double om = 0.1; om < 2; om += 0.1){
        matrix = init_all(m);
        x = iter(matrix, om, eps); 
        arr[j] = it;
        ome[j] = om;
        j++;
        free_m(matrix);
    }
    int min = arr[0];
    double mino = ome[0];
    for(int t = 1; t < 20; t++) {
        if (arr[t] < min) {
            min = arr[t];
            mino = ome[t];
        }
    }
    printf("minimum iteration :%d W: %.1lf\n", min, mino);
    printf("\n");
    return x; 
}
//вычисление нормы дли функции числа обусловленности
double calc_norm(double **m) {
    double max;
    for (int q = 0; q < n; q++) {
        double num = 0;
        for (int p = 0; p < n; p++) {
            num += fabs(m[p][q]);
        }

        if (q == 0 || max <  num){
            max = num;
        }
    }
    return max;
}
//число обусловленности
double condition_number(double **m) {
    printf("CALCULATIN COND NUM\n");
    double **matrix = init_all(m);
    double **ed = init_all(m);
    
    double n1 = calc_norm(matrix);
    ed = get_rev(ed);
    double n2 = calc_norm(ed);
    printf("RESULT: \n%lf\n", n1 * n2);
    free_m(ed);
    free_m(matrix);
    return n1 * n2;
}   

int main(int argc, char const *argv[])
{
    int mgen;
    FILE *f;
    double x;
    sscanf(argv[1], "%d", &mgen);
    if (mgen == 1) {
        n = 4;
        scanf("%lf", &x);
    } else {
        scanf("%d", &n);
        f = fopen("in.txt", "r");
    }
    
    double **matrix = init_m();
    if (mgen) {
        gen_m(matrix, x);
        
    } else {
        get_m(matrix, f);
    }

    double **m_o = orig_g(matrix);
    double **u_o = upgr_g(matrix);
    double det = get_determinant(matrix);
    double **rev = get_rev(matrix);
    double eps;
    scanf("Enter accuracy: %lf", eps);
    double *vec = it_meth(matrix, eps);
    print_v(vec, "result:");

    double c_n = condition_number(matrix);
    free_m(matrix);
    free_m(m_o);
    free(vec);
    free_m(u_o);
    free_m(rev);
    return 0;
}