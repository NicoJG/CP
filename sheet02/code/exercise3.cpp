#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    // read in the a vectors
    vector<Eigen::Vector4d> a = {};
    ifstream file;
    file.open("exercise3.csv");
    if(file.is_open())
    {
        string line;
        while (getline(file,line))
        {
            istringstream iss(line);
            string field;
            Eigen::Vector4d temp;
            int i = 0;
            while (getline(iss,field,','))
            {
                temp[i] = stod(field);
                i++;
            }
            a.push_back(temp);
        }
        
    }

    // a vectors to A matrix
    Eigen::MatrixXd A(4,a.size());
    for(int i=0; i<a.size(); i++)
    {
        A.col(i) = a[i];
    }
    cout << "A:" << endl << A << endl;

    // perform SVD
    Eigen::BDCSVD<Eigen::MatrixXd> svd = A.bdcSvd(Eigen::ComputeFullU);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::VectorXd W = svd.singularValues();

    cout << "U (ONB of the a vetors):" << endl << U << endl;
    cout << "W:" << endl << W << endl;

    return 0;
}