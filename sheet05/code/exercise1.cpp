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
    return sin(x) / x;
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
int main()
{    
    
    std::cout <<"I_1:trap "<< trapezoidal(I_1,1,-1,15) << std::endl;//no idea why b=-1 and a=1 works;N=15 randomly worked...
    std::cout <<"I_1:simps "<< simpson(I_1,-1,1,15) << std::endl;
    std::cout <<"I_2:trap "<< trapezoidal(I_2,0,100,10000000) << std::endl;// error should not be to high with x_max=100 because e^(-100) is already 10^-44
    std::cout <<"I_2:simps "<< simpson(I_2,0,100,10000000) << std::endl;
    std::cout <<"I_3:trap "<< trapezoidal(I_3,-1e8,1e8,15) << std::endl;// I don't know what values for a and b would be apporpriate
    std::cout <<"I_3:simps "<< simpson(I_3,-1e8,1e8,12) << std::endl;
    return 0;   
}