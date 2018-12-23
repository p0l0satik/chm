#include <stdlib.h>
#include <stdio.h>
double func(double x, double y) {

    return -y - (x * x);
}

void runge_kutt2(int steps, double start_y, double start_x, double *sol, double h){
    double y = start_y;

    for(int t = 0; t < steps; t++) {
        double y_wave = y + func(start_x + h * t, y) * h;
        y = y + (func(start_x + h * t, y) + func(start_x + h * (t + 1), y_wave)) * h;
        sol[t + 1] = y; 
    }
}
int main(int argc, char const *argv[])
{
    int steps;
    
    double a, b, ya;
    scanf("enter a, b, ya and steps: %lf%lf%lf%d", &a, &b, &ya, &steps);
    double *sol = calloc(steps, sizeof(sol)); 
    double h = (a - b) / steps;
    runge_kutt2(steps, ya, a, sol, h);
    for (int t = 0; t < steps; t++) {
        printf("%lf %lf\n", a + h * t, sol[t]);
    }
    return 0;
}