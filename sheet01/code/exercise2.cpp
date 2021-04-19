#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    // read in the data
    fstream file("exercise2.csv", ios::in);
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

    

    return 0;
}