#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <fstream>

//Period for harmonic oscillator: T = 2 pi / omega
double const omega = 1.0;
double const T = 2 * M_PI / omega;

//arbitrary vector field
Eigen::VectorXd F(Eigen::VectorXd r)
{
    return -r * omega * omega;
}

//generating y vector for runge kutta
Eigen::VectorXd f(Eigen::VectorXd y)
{
    int n = y.size() / 2;

    Eigen::VectorXd y_n = Eigen::VectorXd(2 * n);

    //from lecture we know : y_vec = (r_vec, v_vec)^T, so each iteration step r_dot and v_dot getting solved (v_vec_dot -> v_vec), (r_vec_dot -> r_vec), therefore each step y has to be reassigned as:
    y_n.head(n) = y.tail(n);
    y_n.tail(n) = F(y.head(n));
    return y_n;
}

//runge kutta 4th order
Eigen::VectorXd rk4(Eigen::VectorXd y, double h)
{
    int n = y.size();

    Eigen::MatrixXd k(n, 4);

    k.col(0) = h * f(y);
    k.col(1) = h * f(y + k.col(0) / 2);
    k.col(2) = h * f(y + k.col(1) / 2);
    k.col(3) = h * f(y + k.col(2));

    y = 1. / 6. * (k.col(0) + 2 * k.col(1) + 2 * k.col(2) + k.col(3));

    return y;
}

//function for solving arbitrary vector field F with rk4
Eigen::MatrixXd ode_solver(Eigen::VectorXd r_0, Eigen::VectorXd v_0, double h, double T)
{
    int N = int(T / h);
    double t = 0;
    int n = r_0.size();

    //saving everything into matrix, after calc matrix can be saved into .log file
    Eigen::MatrixXd outmat(N + 1, int(2 * n + 1));
    Eigen::VectorXd y = Eigen::VectorXd(int(2 * n));

    y << r_0, v_0;

    outmat.row(0) << t, y.transpose();

    for (int i = 0; i < N; i++)
    {
        t += h;
        y += rk4(y, h);
        outmat.row(i + 1) << t, y.transpose();
    }
    return outmat;
}

int main()
{
    std::ofstream outfile;

    //----- Exercise 2.a.1
    //initializing r_0=random and v_0=0
    Eigen::Vector3d r_0 = Eigen::Vector3d::Random(3);
    Eigen::Vector3d v_0 = Eigen::Vector3d::Zero(3);

    //h hasn't to fulfill any condition, so it's randomly chosen.
    double h = 0.1;
    Eigen::MatrixXd outmat = ode_solver(r_0, v_0, h, T);

    outfile.open("build/ex02_a_1.log", std::ios::out);
        outfile << "#t \t x \t y \t z \t vx \t vy \t vz \n";
        outfile << outmat;
    outfile.close();

    //----- Exercise 2.a.1
    //initializing r_0=random and v_0=random and assuming, that Eigen::Random doesn't create parallel vectors
    r_0 = Eigen::Vector3d::Random(3);
    v_0 = Eigen::Vector3d::Random(3);

    outmat = ode_solver(r_0, v_0, h, T);

    outfile.open("build/ex02_a_2.log", std::ios::out);
        outfile << "#t \t x \t y \t z \t vx \t vy \t vz \n";
        outfile << outmat;
    outfile.close();

    //----- Exercise 2.b
    // starting with h=0.1, number of periods=p=10
    h = 0.1;
    double eps = 1e-5;
    double eps_tmp = 1.;
    int p = 10;

    //calc eps_tmp each iteration step, as long as eps_tmp > eps
    while (eps_tmp > eps)
    {
        outmat = ode_solver(r_0, v_0, h, T);
        Eigen::VectorXd r_tmp = outmat.row(outmat.rows() - 1).segment(1, r_0.size());
        eps_tmp = (r_0 - r_tmp).norm();
        h /= 2;
    }
    outfile.open("build/ex02_b.log", std::ios::out);
        outfile << "#eps \t h \n";
        outfile << eps_tmp << "\t" << h;
    outfile.close();

    //----- Exercise 2.c
    //Testing energy conservation for t = 20 via plotting E = 1/2 * v**2 + 1/2 m*w*x**2
    r_0 << 1, 1, 1;
    v_0 << 0, 0, 0;
    h = 0.1;

    double t = 20;

    outmat = ode_solver(r_0, v_0, h, t);

    outfile.open("build/ex02_c.log", std::ios::out);
        outfile << "#t \t x \t y \t z \t vx \t vy \t vz \n";
        outfile << outmat;
    outfile.close();

    return 0;
}