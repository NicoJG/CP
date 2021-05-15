#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>

//setting del_x for derivations
const double f_del_x = 1e-1;

double f(double x)
{
    return x*x - 2;
}

//first and second derivation of function f
double f_primes(double(* func)(double), double x, int i)
{
    double y;
    if (i == 1)
    {
        y = (func(x + f_del_x) - f(x - f_del_x)) / (2 * f_del_x);
    } 
    else if (i == 2)
    {
        y = (func(x + f_del_x) - 2*func(x) + func(x - f_del_x)) / (f_del_x * f_del_x);
    }
    else 
    {
        y = x;
    }

    return y;
}

//bisect method, with simple interval bisection: i.e. d = (a+b)/2
std::tuple<double, int> bisect(double(* func)(double), double a, double b, double c, double eps)
{
    // double min;
    int iter_count = 0;
    while (fabs(c-a) > eps)
    {
        if (c-b > b-a)
        {
            double d = (c + b) / 2.;
            if (func(d) < func(b))
            {   
                a = b; 
                b = d;
            }
            else
            {
                c = d;
            }
        }
        else
        {
            double d = (b + a) / 2.;
            if (func(d) < func(b))
            {
                a = d;
            }
            else
            {
                c = b;
                b = d;
            }
        }
        iter_count++;
    }

    double min = (c-a)/2;
    return std::make_tuple(min, iter_count);
}

//newton method
std::tuple<double, int> newton(double(* func)(double), double x_0, double eps)
{
    double del_x = -f_primes(func, x_0, 1) / f_primes(func, x_0, 2);
    int iter_count = 0;

    while (fabs(del_x) > eps)
    {
        x_0 += del_x;
        del_x = -f_primes(func, x_0, 1) / f_primes(func, x_0, 2);
        iter_count++;
    }

    return std::make_tuple(x_0, iter_count);
}

int main()
{
    double a, b, c, x_0, eps;
    a = -0.5;
    b = -0.1;
    c = 2.;
    x_0 = 1.;
    eps = 1e-9;

    double min_bisect, min_newton;
    int iter_bisect, iter_newton;
    
    std::tie(min_bisect, iter_bisect) = bisect(f, a, b, c, eps);
    std::tie(min_newton, iter_newton) = newton(f, x_0, eps);

    std::ofstream outfile;

    outfile.open("build/ex03.log", std::fstream::out);
        outfile <<  "#i_bisect    min_bisect    i_newton    min_newton" << "\n";
        outfile << iter_bisect << "\t" << min_bisect << "\t" << iter_newton << "\t" << min_newton;
    outfile.close();

    return 0;
}