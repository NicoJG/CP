#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"

#include "profiler.h"


using namespace std;

int main()
{   //a)
    int N = 3;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(N, N);
    Eigen::VectorXd b = Eigen::VectorXd::Random(N);
    cout << "Here is the matrix m:" << endl << m << endl;
    cout << "Here is the vector b:" << endl << b << endl;

    Eigen::FullPivLU<Eigen::MatrixXd> PPLU(m);
    Eigen::VectorXd x = PPLU.solve(b);
    if ((m * x).isApprox(b))
    {
        cout << "Here is a solution x to the equation mx=b:" << endl << x << endl;
    }
    else
        cout << "The equation mx=y does not have any solution." << endl;
    //b)
    Profiler::init(3);// initializing 3 profiler
    vector<double> t_step1(1000);
    vector<double> t_step2(1000);
    vector<double> t_step3(1000);
    int i = 1;
    for (i; i < 1000; i++)
    {
        Profiler::resetAll();// reset all profiler at beginning of loop
        //step 1
        Profiler::start(0);
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(i, i);
        Eigen::VectorXd y = Eigen::VectorXd::Random(i);
        Profiler::stop(0);
        //step 2
        Profiler::start(1);
        Eigen::FullPivLU<Eigen::MatrixXd> LU(A);
        Profiler::stop(1);
        //step 3
        Profiler::start(2);
        Eigen::VectorXd t = LU.solve(y);
        Profiler::stop(2);
        // saving times in the 3 vectors
        t_step1[i] = Profiler::getTimeInS(0);
        t_step2[i] = Profiler::getTimeInS(1);
        t_step3[i] = Profiler::getTimeInS(2);
    }

    auto file = fopen("times.dat", "w");
    fprintf(file, "# t1 t2 t3\n");
    for (int q = 1; q < 1000; q++) {
        fprintf(file, "%.10f %.10f %.10f\n", t_step1[q], t_step2[q], t_step3[q]);

    }
    fclose(file);
}
