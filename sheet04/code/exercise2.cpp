#include <iostream>
#include <math.h>
#include <fstream> 

//potential with first charge distribution
double phi_1(int n, double x_p, double y_p, double z_p)
{
    return 1 / sqrt(pow(0.1*n - x_p, 2) + y_p*y_p + z_p*z_p + 1e-6);
}

//potential with second charge distribution
double phi_2(int n, double x_p, double y_p, double z_p)
{   
    return x_p / sqrt((0.1*n-x_p)*(0.1*n-x_p) + y_p*y_p + z_p*z_p + 1e-6);
}

//midpoint rule for calculating the integral
double midpoint(double (* func)(int, double, double, double), double a, double b, int N, int n)
{
    double h = (b - a) / N;
    double result = 0;
    
    for (int i = 1; i < N; i++)
    {
            for (int j = 1; j < N; j++)
        {
                for (int k = 1; k < N; k++)
            {
                result += func(n, a + i*h - (h/2), a + j*h - (h/2), a + k*h - (h/2));
            }
        }
        
    }
    return h*h*h*result;
}


int main()
{
    double a = -1;
    double b = 1;
    int N = 100;

    std::ofstream outfile1;
    std::ofstream outfile2;

    outfile1.open("build/ex02_a.log", std::fstream::out);
    outfile2.open("build/ex02_b.log", std::fstream::out);

    outfile1 << "# n" << "\t" << "phi_1" <<"\n"; 
    outfile2 << "# n" << "\t" << "phi_2" <<"\n"; 

    //calculating and saving all values (inside and outside the cube)
    for (int n = 0; n < 81; n++)
        {   
            outfile1 << 0.1*n << "\t" <<  midpoint(phi_1, a, b, N, n) << "\n";
            outfile2 << 0.1*n << "\t" <<  midpoint(phi_2, a, b, N, n) << "\n";
        }
    outfile1.close();
    outfile2.close();
    return 0;
}