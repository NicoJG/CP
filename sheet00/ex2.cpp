#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

int main() 
{
    // a)
    float a_min = 1000;
    float a_step = 100;
    float a_max = 1000000;

    ofstream file_a;
    file_a.open("build/ex2_a.csv");
    file_a << "x,y_exact,y_rounded,y_better,rel_err_rounded,rel_err_better" << endl;

    for(float x=a_min;x<=a_max;x+=a_step)
    {
        double y_exact = 1./sqrt((double)x)-1/sqrt((double)x+1.);
        float y_rounded = 1./sqrt(x)-1/sqrt(x+1.);
        float y_better = (sqrt(x+1)-sqrt(x)) / (sqrt(x)*sqrt(x+1));
        double rel_err_rounded = ((double)y_rounded - y_exact) / y_exact;
        double rel_err_better = ((double)y_better - y_exact) / y_exact;

        file_a << x << "," << y_exact << "," << y_rounded << "," << y_better << "," << rel_err_rounded << "," << rel_err_better << endl;
    } 

    file_a.close();

    // b)
    
    float b_min = 0.000001;
    float b_step = 0.00001;
    float b_max = 0.002;

    ofstream file_b;
    file_b.open("build/ex2_b.csv");
    file_b << "x,y_exact,y_rounded,y_better,rel_err_rounded,rel_err_better" << endl;

    for(float x=b_min;x<=b_max;x+=b_step)
    {
        double y_exact = (1.-cos((double)x))/sin((double)x);
        float y_rounded = (1.-cos(x))/sin(x);
        float y_better = tan(x/2.);
        double rel_err_rounded = ((double)y_rounded - y_exact) / y_exact;
        double rel_err_better = ((double)y_better - y_exact) / y_exact;

        file_b << x << "," << y_exact << "," << y_rounded << "," << y_better << "," << rel_err_rounded << "," << rel_err_better << endl;
    } 

    file_b.close();
    
    // c)
    
    float c_x_min = 1;
    float c_x_step = 10000;
    float c_x_max = 1000000;
    float c_delta_min = 0.000001;
    float c_delta_step = 0.0001;
    float c_delta_max = 0.02;

    ofstream file_c;
    file_c.open("build/ex2_c.csv");
    file_c << "x,delta,y_exact,y_rounded,y_better,rel_err_rounded,rel_err_better" << endl;

    for(float x=c_x_min;x<=c_x_max;x+=c_x_step)
    {
        for(float delta=c_delta_min;delta<=c_delta_max;delta+=c_delta_step)
        {
            double y_exact = sin((double)x+(double)delta)-sin((double)x);
            float y_rounded = sin(x+delta)-sin(x);
            float y_better = sin(x+delta)-sin(x);
            double rel_err_rounded = ((double)y_rounded - y_exact) / y_exact;
            double rel_err_better = ((double)y_better - y_exact) / y_exact;

            file_c << x << "," << delta << "," << y_exact << "," << y_rounded << "," << y_better << "," << rel_err_rounded << "," << rel_err_better << endl;
        }
    } 

    file_c.close();
    
    return 0;
}
