#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <vector>

using namespace std;
//setting h for derivations
const double h = 1e-1;
double rosenbrock(double x, double y, double lambda)
{
    return pow(1-x*lambda,2.0) + 100. * pow(y*lambda - pow(x*lambda,2.), 2.0);
}

double f_b(double x, double y, double lambda)
{
    return pow(1.+(exp(-10.* pow(x*y*lambda*lambda -3.,2.))/(pow(x*lambda,2.)+pow(y*lambda,2.))),-1.);
}
//first and second partial derivation of function f for lambda 
double f_primes(double(* func)(double, double , double), double x, double y, double lambda,  int i)
{
    double r;
    if (i == 1)
    {
        r = (func(x, y, lambda + h) - func(x, y, lambda - h)) / (2 * h);
    } 
    else if (i == 2)
   {
        r = (func(x, y, lambda + h) - 2*func(x, y, lambda) + func(x, y, lambda - h)) / (h * h);
    }
    else 
    {
        r = x;
    }

    return r;
}

//no stopping criteria just 10000 steps because the algorithm took too long otherwise
double newton(double(* func)(double, double, double), double x_0, double y_0, double lambda)
{
    double del_l = -f_primes(func, x_0, y_0, lambda, 1) / f_primes(func, x_0, y_0, lambda, 2);
    //int iter_count = 0;
    //double eps = 1e-3;
    //while (fabs(del_l) > eps)
    for (int i = 0; i < 10000; i++)
    {
        lambda += del_l;
        del_l = -f_primes(func, x_0, y_0, lambda, 1) / f_primes(func, x_0, y_0, lambda, 2);
        //iter_count++;
    }

    return lambda;
}


//gradient descent
std::tuple<vector<double>, vector<double>> grad(double(* func)(double, double, double),double x_0, double y_0, double lambda, int iter_count)//, double g_c)
{   double g_x=1;
    double g_y=1;
    vector<double> x_n;
    x_n.push_back( x_0 );
    vector<double> y_n;
    y_n.push_back( y_0 );
    int i = 0;

    //while(sqrt(g_x*g_x + g_y*g_y)> g_c)
    for (int i = 0; i < iter_count; i++)
    {
    
    g_x = - (func(x_n[i] + h, y_n[i], 1.) - func(x_n[i] - h, y_n[i], 1.))/2*h;//calculating the gradient elemntwise
    g_y = - (func(x_n[i], y_n[i] + h, 1.) - func(x_n[i], y_n[i] - h, 1.))/2*h;
    lambda = newton(func,x_n[i]/lambda + g_x , y_n[i]/lambda + g_y, lambda);//calculating stepsize;cant find segmentation error 
    
    x_n.push_back(x_n[i] + lambda * g_x);//saving new x_i & y_i
    y_n.push_back(y_n[i] + lambda * g_x);
    
    }
    return std::make_tuple(x_n, y_n);
}
//Conjugate gradient method
std::tuple<vector<double>, vector<double>> cgrad(double(* func)(double, double, double),double x_0, double y_0, double lambda, int iter_count)//, double g_c)
{
    vector<double> x_n;
    x_n.push_back( x_0 );
    vector<double> y_n;
    y_n.push_back( y_0 );
    
    double g_x = - (func(x_n[0] + h, y_n[0], 1.) - func(x_n[0] - h, y_n[0], 1.))/2*h;//calculate g_0
    double g_y = - (func(x_n[0], y_n[0] + h, 1.) - func(x_n[0], y_n[0] - h, 1.))/2*h;

    double p_x = g_x;
    double p_y = g_y;
    //double lambda = 0.5; 

    for (int i = 0; i < iter_count; i++)
    {   
        lambda = newton(func,x_n[i]/lambda + g_x , y_n[i]/lambda + g_y, lambda);
        x_n.push_back(x_n[i]+ lambda * p_x);//calculate x_i and store it in x_n
        y_n.push_back(y_n[i]+ lambda * p_y);

        double g_xprior = g_x;
        double g_yprior = g_y;
        
        g_x = - (func(x_n[i+1] + h, y_n[i+1], 1.) - func(x_n[i+1] - h, y_n[i+1],1.))/2*h;//calculate g_i+1
        g_y = - (func(x_n[i+1], y_n[i+1] + h, 1.) - func(x_n[i+1], y_n[i+1]+h , 1.))/2*h;
        
        double mu = (g_x*g_x + g_y*g_y)/(g_xprior*g_xprior + g_yprior*g_yprior);//calculate new minimization direction
        p_x = g_x + mu * p_x;
        p_y = g_y + mu * p_y;
        
    }
    return std::make_tuple(x_n, y_n);
}
int main()
{   
    //double g_c = 1e-2;
    vector<double> x_g, y_g, x_cg, y_cg;

    std::tie(x_g, y_g) = grad(rosenbrock,-1.,-1., 1., 1000);//,g_c);
    
    auto file1 = fopen("build/ex02gd.dat", "w");
    fprintf(file1, "# x_n y_n\n");
    for(int i = 0; i < x_g.size(); i++){
        fprintf(file1, "%.10f %.10f\n", x_g[i], y_g[i]);
    }
    fclose(file1);

    std::tie(x_cg, y_cg) = cgrad(rosenbrock,-1.,-1., 1.,7);
    
    auto file2 = fopen("build/ex02cgd.dat", "w");
    fprintf(file2, "# x_n y_n\n");
    for(int i = 0; i < x_cg.size(); i++){
        fprintf(file2, "%.10f %.10f\n", x_cg[i], y_cg[i]);
    }
    fclose(file2);
    
    vector<double> x_1, y_1, x_2, y_2, x_3, y_3;
    
    std::tie(x_1, y_1) = cgrad(rosenbrock,1.5,2.3, 1.,5);
    cout<<"X_min 1 :"<< endl << x_1.back() << endl << y_1.back() << endl;
    std::tie(x_2, y_2) = cgrad(rosenbrock,-1.7,-1.9, 1.,20);
    cout<<"X_min 2 :"<< endl << x_2.back() << endl << y_2.back() << endl;
    std::tie(x_3, y_3) = cgrad(rosenbrock,0.5,0.6, 1.,20);
    cout<<"X_min 3 :"<< endl << x_3.back() << endl << y_3.back() << endl;
    return 0;
}