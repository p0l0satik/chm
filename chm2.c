#include <stdlib.h>
#include <stdio.h>
#include <math.h>
double func(double x, double y) {
    return -y - (x * x);
}

double sys(double x, double v, double u, int f) {
    if (f == 0) {
        return u + x -2;
    } 
    if (f == 1) {
        return  5 * u + 6;
    }
    printf("wrong input\n");
    return -1;
}

double sys1(double x, double v, double u, int f) {
    if (f == 0) {
        return sin(x + u) - 1.1 * v;
    } 
    if (f == 1) {
        return 2.5 * u - (x + v) * (x + v);
    }
    printf("wrong input\n");
    return -1;
}


void get_func(int steps, double start_x, double *sol, double h) {
    double x =  start_x;
    for (int t = 0; t < steps; t++) {
        sol[t] = -x * x + 2 * x - 2 + 12 * exp(-x);
        x += h;
    }
}

void runge_kutt2(int steps, double start_y, double start_x, double *sol, double h){
    double y = start_y;
    double x = start_x;
    sol[0] = y;
    for(int t = 0; t < steps; t++) {
        double y_wave = y + func(x, y) * h;
        y = y + (func(x, y) + func(x + h, y_wave)) * h / 2;
        sol[t + 1] = y; 
        x += h;
    }
}

void runge_kutt_sys2(int steps, double *y, double x, double **sol, double h) {
    sol[0][0] = y[0];
    sol[1][0] = y[1];
    for (int t = 0; t < steps; t++){
        double y_w[2];
        for (int k = 0; k < 2; k++){
            y_w[k] = y[k] + sys(x, y[0], y[1], k) * h;
        }
        for (int k = 0; k < 2; k++){
            if (k == 0) {
                y[0] = y[0] + (sys(x, y[0], y[1], k) + sys(x + h, y_w[0], y[1], k)) * h / 2;
            } else {
                y[1] = y[1] + (sys(x, y[0], y[1], k) + sys(x + h, y[0], y_w[1], k)) * h / 2;
            }
            
            sol[k][t + 1] = y[k];
        }
        x += h;
    }
}

void runge_kutt_sys4(int steps, double *y, double x, double **sol, double h) {
    for (int t = 0; t < steps; t++){
        for (int k = 0; k < 2; k++){
            double k1 = sys(x, y[0], y[1], k);
            double k2 = sys(x + h / 2, y[0] + h / 2 * k1, y[1] + h / 2 * k1, k);
            double k3 = sys(x + h / 2, y[0] + h / 2 * k2, y[1] + h / 2 * k2, k);
            double k4 = sys(x + h, y[0] + h * k3, y[1] + h * k3, k);
            y[k] = y[k] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
            sol[k][t + 1] = y[k]; 
        }
        x += h;
    }
}

void runge_kutt4(int steps, double start_y, double start_x, double *sol, double h){
    double y = start_y;
    double x = start_x;
    sol[0] = y;
    for(int t = 0; t < steps; t++) {
        double k1 = func(x, y);
        double k2 = func(x + h / 2, y + h / 2 * k1);
        double k3 = func(x + h / 2, y + h / 2 * k2);
        double k4 = func(x + h, y + h * k3);
        y = y + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
        sol[t + 1] = y; 
        x += h;
    }
}
int main(int argc, char const *argv[])
{
    int steps;
    
    double a, b;
    printf("enter a, b,and steps: \n");
    scanf("%lf%lf%d", &a, &b, &steps);
    double h = (b - a) / steps;
    if (argv[0][0] == '1') {
        double ya;
        printf("enter y\n");
        scanf("%lf", &ya);
        double *sol = calloc(steps, sizeof(sol)); 
        double *orig = calloc(steps, sizeof(sol));
        runge_kutt2(steps, ya, a, sol, h);
        runge_kutt4(steps, ya, a, sol, h);
        get_func(steps, a, orig, h);
        FILE *twout = fopen("2out.txt", "w");
        FILE *frout = fopen("4out.txt", "w");
        for (int t = 0; t < steps; t++) {
            fprintf(twout, "(%lf;%lf) ", a + h * t, sol[t]);
        }
        for (int t = 0; t < steps; t++) {
            fprintf(frout, "(%lf;%lf) ", a + h * t, sol[t]);
        }
    }

    double **sys_sol = calloc(2, sizeof(sys_sol));
    for (int t = 0; t < 2; t++) {
        sys_sol[t] = calloc(steps, sizeof(sys_sol[t]));
    }
    double y[2];
    printf("enter y1 y2\n");
    scanf("%lf%lf", &y[0], &y[1]);
    runge_kutt_sys2(steps, y, a, sys_sol, h);
    FILE *out[2] = {fopen("s1out.txt", "w"), fopen("s2out.txt", "w")};
    for (int k = 0; k < 2; k++) {
        for (int t = 0; t < steps; t++) {
            fprintf(out[k], "(%lf;%lf) ", a + h * t, sys_sol[k][t]);
        }
    }
        
    return 0;
}