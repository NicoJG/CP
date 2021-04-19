#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>

using namespace std;

int main()
{
    Eigen::Vector3d a1 = {0.5,sqrt(3.)/2.,0};
    Eigen::Vector3d a2 = {-0.5,sqrt(3.)/2.,0};
    Eigen::Vector3d a3 = {0,0,1};
    Eigen::Matrix3d A;
    A << a1,a2,a3;

    // b)
    Eigen::PartialPivLU<Eigen::Matrix3d> lu = A.lu();
    Eigen::Matrix3d P = lu.permutationP();
    Eigen::Matrix3d LU = lu.matrixLU();

    Eigen::Vector3d b = {2,0,2};
    Eigen::Vector3d x = lu.solve(b);
    
    cout << "A:" << endl << A << endl;
    cout << "P:" << endl << P << endl;
    cout << "LU:" << endl << LU << endl;
    cout << "x:" << endl << x << endl;

    // c)
    Eigen::Vector3d c = {1,2.*sqrt(3.),3};
    Eigen::Vector3d y = lu.solve(c);

    cout << "y:" << endl << y << endl;

    // d)
    Eigen::Matrix3d A2;
    A2 << a3,a2,a1;

    Eigen::PartialPivLU<Eigen::Matrix3d> lu2 = A2.lu();
    Eigen::Matrix3d P2 = lu2.permutationP();
    Eigen::Matrix3d LU2 = lu2.matrixLU();

    cout << "A2:" << endl << A2 << endl;
    cout << "P2:" << endl << P2 << endl;
    cout << "LU2:" << endl << LU2 << endl;

    return 0;
}