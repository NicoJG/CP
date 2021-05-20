#include <iostream>
#include <math.h>

double I_1(double x)
{
    return exp(-x) / x;
}
double I_2(double x)
{   if(x==0){
    return 0;
}   else{
    return exp(-x) / (sqrt(x));
}
    
}
double I_3(double x)
{
    return sin( x ) / x;
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

int main()
{    
    
    std::cout <<"I_1:trap "<< trapezoidal(I_1,1,-1,15) << std::endl;
    
    std::cout <<"I_2:trap "<< trapezoidal(I_2,0,100,1e7) << std::endl;// error should not be to high with x_max=100 because e^(-100) is already 10^-44
    
    std::cout <<"I_3:mid "<< 2 * (midpoint(I_3,0,1,1e5) + midpoint(I_3,1,10000,1e5)); //midpoint bc its belongs to the OCN and ignores the singularity
    
    return 0;   
}