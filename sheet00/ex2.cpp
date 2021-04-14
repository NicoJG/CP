#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

float a(float x)
{
    return 1 / sqrt(x) - 1 / sqrt(x + 1);
}

int main() 
{
    // a)
    vector<float> x0 = {};
    for(float i=1000;i<1000000;i+=100)
    {
        x0.push_back(i);
    }

    ofstream file0;
    file0.open("build/ergebnisse1.csv");

    vector<float> y0 = {};
    
    for(float y:x0)
    {
        file0 <<  y << "," << a(y) << endl;
        y0.push_back(a(y));
    }
    file0.close();

        vector<double> x1 = {};
    for(double i=1000;i<1000000;i+=100)
    {
        x1.push_back(i);
    }

    ofstream file1;
    file1.open("build/exakt1.csv");

    vector<float> y1 = {};
    
    for(double y:x1)
    {
        file1 <<  y << "," << a(y) << endl;
        y1.push_back(a(y));
    }
    file1.close();

    ofstream diff;
    diff.open("build/diff.csv");

    for(int i=0; i <= y0.size(); i++)
    {
        diff << y0[i]-y1[i] << endl;
    }
    
    diff.close();

    


    // b)

    
    return 0;
}
