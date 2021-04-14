#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

float a0(float x)
{
    return (1-cos(x))/sin(x);
}

double a1(double x)
{
    return (1-cos(x))/sin(x);
}


int main() 
{
    // a)

    //float

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
        file0 <<  y << "," << a0(y) << endl;
        y0.push_back(a0(y));
    }
    file0.close();

    //double

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
        file1 <<  y << "," << a1(y) << endl;
        y1.push_back(a1(y));
    }
    file1.close();


    //difference

    ofstream diff;
    diff.open("build/diff.csv");

    for(int i=0; i <= y0.size(); i++)
    {
        diff << (double)y0[i]-y1[i]/(y1[i]) << endl;
    }
    
    diff.close();

    


    // b)

    
    return 0;
}
