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
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    A << x.dot(x), x.sum(), 
        x.sum(), x.size();
    b << x.dot(y), y.sum();

    Eigen::PartialPivLU<Eigen::Matrix2d> lu = A.lu();

    Eigen::Matrix2d P = lu.permutationP();
    Eigen::Matrix2d LU = lu.matrixLU();

    Eigen::Vector2d mn = lu.solve(b);

    cout << "A:" << endl << A << endl;
    cout << "P:" << endl << P << endl;
    cout << "LU:" << endl << LU << endl;
    cout << "b:" << endl << b << endl;
    cout << "mn:" << endl << mn << endl;

    // save m and n for Python
    ofstream file2("build/exercise2_mn.csv");
    if(file2.is_open())
    {
        file2 << "m,n" << endl;
        file2 << mn(0) << "," << mn(1) << endl;
    }
    file2.close();

    return 0;
}