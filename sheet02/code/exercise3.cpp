#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <string>
#include "read_write_csv.cpp"

using namespace std;

int main()
{
    // read in the A matrix
    Eigen::MatrixXd A = readCSV("exercise3.csv", false).transpose();
    cout << "A:" << endl << A << endl;
    writeCSV("build/exercise3_A.csv",A,0);

    // perform SVD
    Eigen::BDCSVD<Eigen::MatrixXd> svd = A.bdcSvd(Eigen::ComputeFullU);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::VectorXd W = svd.singularValues();

    cout << "U (ONB of the a vetors):" << endl << U << endl;
    writeCSV("build/exercise3_U.csv",U,5);
    cout << "W:" << endl << W << endl;
    writeCSV("build/exercise3_W.csv",W,5);

    return 0;
}