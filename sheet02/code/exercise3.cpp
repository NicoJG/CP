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
    return 0;
}