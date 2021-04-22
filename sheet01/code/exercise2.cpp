#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    // read in the data
    ifstream file("exercise2.csv", ios::in);
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    if(file.is_open())
    {
        string line;
        int i = 0;
        getline(file,line);
        while (getline(file,line))
        {
            istringstream iss(line);
            string field;
            while(getline(iss,field,','))
            {
                x.conservativeResize(i+1);
                x[i] = stod(field);
                getline(iss,field,',');
                y.conservativeResize(i+1);
                y[i] = stod(field);
                i++;
                
            }
        }
    }
    file.close();

    cout << "x:" << endl << x << endl;
    cout << "y:" << endl << y << endl;
    
    // c)
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    Eigen::MatrixXd A(x.size(),2);
    A << x, ones;
    Eigen::MatrixXd A_T = A.transpose();
    Eigen::MatrixXd P = A_T*A;
    Eigen::VectorXd b = A_T * y;

    cout << "A:" << endl << A << endl;
    cout << "P:" << endl << P << endl;
    cout << "b:" << endl << b << endl;

    Eigen::PartialPivLU<Eigen::MatrixXd> lu = P.lu();

    Eigen::MatrixXd P_LU = lu.permutationP();
    Eigen::MatrixXd LU = lu.matrixLU();

    cout << "P_LU:" << endl << P_LU << endl;
    cout << "LU:" << endl << LU << endl;

    Eigen::VectorXd m = lu.solve(b);

    cout << "m:" << endl << m << endl;

    // save m and n for Python
    ofstream file2("build/exercise2_mn.csv");
    if(file2.is_open())
    {
        file2 << "m,n" << endl;
        file2 << m(0) << "," << m(1) << endl;
    }
    file2.close();

    return 0;
}