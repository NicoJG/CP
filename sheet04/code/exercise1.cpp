#include <iostream>
#include <math.h>


double f(double x)
{
    return exp(-x) / x;
}

double g(double x)
{
    return x * sin(1 / x);
}

double midpoint(double( *func)(double), double a, double b, int N)
{   
    double h = (b - a) / N;
    double result = 0;

    for (int i = 1; i < N; i++)
    {
        result += func(a + i*h - h/2);
    }

    return h * result;
}

double trapezoidal(double(* func)(double), double a, double b, int N)
{
    double h = (b - a) / N;
    double result = func(a)/2 + func(b)/2;
    double x = a + h;

    for (int i = 1; i < N; i++)
    {
        result += func(x);
        x += h;
    }
    return h * result;
}

double simpson(double(* func)(double), double a, double b, int N)
{
    double h = (b - a) / N;
    double result = func(a) + 4*func(a + h/2) + func(b);
    double x = a + h;

    for (int i = 1; i < N; i++)
    {
        result += 2*func(x) + 4*func(x + h/2);
        x += h;
    }
    return h / 6 * result;
}

double integrate(double(* func_int)(double (* func)(double), double, double, int), double(* func)(double), double a, double b,  double err_abs)
{
    int N = 2;
    double result_0 = func_int(func, a, b, N);
    double result_1 = 0;
    double err_tmp = 1;
    int iterator = 0;

    while(err_tmp > err_abs)
    {
        N *= 2;

        result_1 = func_int(func, a, b, N);
        err_tmp = abs(result_1 - result_0)/result_1;
        result_0 = result_1;

        iterator += 1;
    }

    std::cout << "i: " << iterator << ", Ergebnis: " <<result_0 << std::endl;

    return result_0;
}

int main()
{
    double exact_f = 0.219384;
    double exact_g = 0.37853;
    double err = 1e-4;
    //adding err to lower limit, because of isolated singularity
    double I_1_mid = integrate(midpoint, f, 1, 100, err);
   
    double I_2_mid = integrate(midpoint, g, 0+err, 1, err);
    
    double I_1_trap = integrate(trapezoidal, f, 1, 100, err);
    
    double I_2_trap = integrate(trapezoidal, g, 0+err, 1, err);
    
    double I_1_simp = integrate(simpson, f, 1, 100, err);
    
    double I_2_simp = integrate(simpson, g, 0+err, 1, err);
    

    return 0;
}