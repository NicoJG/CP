#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

// a)
// calculates the eigenfrequencies of a linear chain of oscillators (open at the ends)
// with masses m_i, spring constants k_i, equillibrium lengths l_i
// N m_i but only N-1 k_i and l_i
VectorXcd eigenfrequencies(int N, VectorXd m, VectorXd k, VectorXd l)
{
    MatrixXd A = MatrixXd::Zero(N, N);
    // first and last rows are special
    A.row(0).head(2) << -k(0) / m(0), k(0) / m(0);
    A.row(N - 1).tail(2) << k(N - 2) / m(N - 1), -k(N - 2) / m(N - 1);
    // fill the matrix
    for (int i = 1; i < N - 1; i++)
    {
        A.row(i).segment(i - 1, 3) << k(i - 1) / m(i), -(k(i - 1) + k(i)) / m(i), k(i) / m(i);
    }
    VectorXcd lambda = A.eigenvalues();
    cout << "Matrix A" << endl << A << endl;
    return (-lambda).cwiseSqrt();

}

int main()
{
    // b)
    int N = 10;
    VectorXd m = VectorXd::LinSpaced(N, 1, N);
    VectorXd k = VectorXd::LinSpaced(N - 1, N - 1, 1);
    VectorXd l = VectorXd::LinSpaced(N - 1, 1, N - 1);
    l.array() = (5 - l.array()).abs() + 1;

    cout << "m:" << endl << m.transpose() << endl;
    cout << "k:" << endl << k.transpose() << endl;
    cout << "l:" << endl << l.transpose() << endl;

    cout << "eigenfrequencies:" << endl << eigenfrequencies(N, m, k, l).transpose() << endl;
    return 0;
}