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

    // a)
    Eigen::Vector3d b = {2,0,2};

    Eigen::PartialPivLU<Eigen::Matrix3d> lu = A.lu();
    Eigen::Matrix3d P = lu.permutationP();
    Eigen::Matrix3d LU = lu.matrixLU();

    
    Eigen::Vector3d x = lu.solve(b);
    
    cout << "A:" << endl << A << endl;
    cout << "P:" << endl << P << endl;
    cout << "LU:" << endl << LU << endl;
    cout << "x:" << endl << x << endl;

    // b)
    return 0;
}