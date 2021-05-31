#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

//setting h for derivations
const double h = 1e-1;

double f(double x1, double x2) //rosenblut function; giving the vector coeffwise and not as a vector to calculate the gradient in a simpler way
{
    return pow(1 - x1, 2.) + 100 * pow(x2 - pow(x1, 2.), 2.);
}
VectorXd gradient(double (*func)(double, double), double x1, double x2) //calculating the gradient of a function for a certain vector; not dynamic that way; only functions with 2 arguments
{
    VectorXd g(2);
    g(0) = -(func(x1 + h, x2) - func(x1 - h, x2)) / 2 * h;
    g(1) = -(func(x1, x2 + h) - func(x1, x2 - h)) / 2 * h;

    return g;
}
double backlinesearch(double (*func)(double, double), double x1, double x2, VectorXd g) //backtracking linesearch; taken from the german wikipedia search because our implementation of the stepsize calculation from the last sheet did not work
{
    double alpha = 1.;
    double sigma = 1e-4;
    double rho = 0.5; //factor by which alpha is reduced
    int i = 0;
    RowVectorXd t = -1. * g.transpose(); //transposed directional derivative
    while (func(x1 + alpha * g(0), x2 * alpha * g(1)) > func(x1, x2 + sigma * alpha * t * g))
    {
        alpha = rho * alpha;
        i++;
        if (i == 100)
        {
            break; // condition if line search does not converge
        }
    }
    return alpha;
}

std::tuple<VectorXd, int> bfgs(double (*func)(double, double), double x1, double x2, MatrixXd C0, double epps)
{
    VectorXd g = gradient(func, x1, x2);            // calculating b_0 (step 1 in the script)
    double alpha = backlinesearch(func, x1, x2, g); // step 2
    Vector2d x0;
    x0 << x1, x2;
    VectorXd x_new = x0 + alpha * g;
    VectorXd g_new = gradient(func, x_new(0), x_new(1));
    VectorXd S_0 = x_new - x0;
    VectorXd y_0 = g_new - g;
    //step 3 will be done outside of the function(only the first calculation of c0)
    VectorXd p_k;
    VectorXd x_prior;
    VectorXd g_prior;
    int i = 0;
    while (g_new.norm() > epps)
    {
        MatrixXd I = Matrix<double, 2, 2>::Identity();
        double rho_k = 1 / (y_0.transpose() * S_0);
        C0 = (I - rho_k * S_0 * y_0.transpose()) * C0 * (I - rho_k * y_0 * S_0.transpose()) + rho_k * S_0 * S_0.transpose(); //step 4
        p_k = -C0 * g_new;
        x_prior = x_new;
        g_prior = g_new;
        x_new = x_new + p_k;
        g_new = gradient(func, x_new(0), x_new(1));
        S_0 = x_new - x_prior;
        y_0 = g_new - g_prior;
        i++;
    }
    return std::make_tuple(x_new, i);
}

int main()
{
    double x1, x2;
    x1 = x2 = -1.;
    double epps = 1e-5;
    MatrixXd I = Matrix<double, 2, 2>::Identity();
    int iter_C1, iter_C2,iter_C3;
    Vector2d min_c1,min_c2,min_c3;

    //exact inverse of Hesse Matrix
    Matrix2d H1;//hessian matrix of rosenblut function
    Matrix2d C1;//inverse of H1
    H1 << 1200 * pow(x1, 2.) - 400 * x2 + 2, -400 * x1,
        -400 * x1, 200;
    double det = 1 / H1.determinant();
    C1 << 200, 400 * x1,
        400 * x1, 1200 * pow(x1, 2.) - 400 * x2 + 2;//calculated inverse analyticly
    C1 = det * C1;
    std::tie(min_c1, iter_C1) = bfgs(f,x1,x2,C1,epps);
    cout << "min_c1:  "<< endl << min_c1 << endl << "iter_c1:  " << iter_C1<< endl;
    //diagonal matrix on whose diagonal are the inverse diagonal elements
    Matrix2d C2;
    C2 << pow(1200 * pow(x1, 2.) - 400 * x2 + 2, -1.), 0,
        0, pow(200, -1.);
    std::tie(min_c2, iter_C2) = bfgs(f,x1,x2,C2,epps);
    cout << "min_c2:  "<< endl << min_c2 << endl << "iter_c2:  " << iter_C2<< endl;
    //multiple of the unit matrix with a prefactor in the "typical order of magnitude"
    Matrix2d C3;
    C3 = f(x1, x2) * I;
    std::tie(min_c3, iter_C3) = bfgs(f,x1,x2,C3,epps);
    cout << "min_c3:  "<< endl << min_c3 << endl << "iter_c3:  " << iter_C3<< endl;
    return 0;
}