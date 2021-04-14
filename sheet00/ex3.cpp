#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;


const double delta_t = 0.01;
const int max_n = 10/delta_t;

vector<vector<double>> euler(double anfang_norm, double anfang_sym_0, double anfang_sym_1)
{
    double x = anfang_norm;
    vector<double> y_norm = {};
    for( int n=0 ; n <= max_n ; n++ )
    {
        x = x * ( 1 - delta_t);
        y_norm.push_back(x);
    }
    
    
    vector<double> y_sym = {anfang_sym_0,anfang_sym_1};
    for(int n=1; n < max_n; n++)
    {
        y_sym.push_back(-2 * delta_t * y_sym[n] + y_sym[n-1]);
    }

    vector<vector<double>> result = {};
    result.push_back(y_norm);
    result.push_back(y_sym);
    return result;
}

int main()
{
    vector<vector<double>> result_a = euler(1,1,exp(-delta_t));
    vector<double> y_norm_a = result_a[0];
    vector<double> y_sym_a = result_a[1];

    vector<vector<double>> result_b = euler(1-delta_t,1,1-delta_t);
    vector<double> y_norm_b = result_b[0];
    vector<double> y_sym_b = result_b[1];
    
    vector<double> y_exact = {};
    for(int n=0; n <= max_n; n++)
    {
        y_exact.push_back(exp(-n*delta_t));
    }



    vector<double> rel_err_norm_a = {};
    vector<double> rel_err_sym_a = {};
    for(int n=0; n <= max_n; n++)
    {
        rel_err_norm_a.push_back((y_norm_a[n]-y_exact[n])/y_exact[n]);
        rel_err_sym_a.push_back((y_sym_a[n]-y_exact[n])/y_exact[n]);
    }

    vector<double> rel_err_norm_b = {};
    vector<double> rel_err_sym_b = {};
    for(int n=0; n <= max_n; n++)
    {
        rel_err_norm_b.push_back((y_norm_b[n]-y_exact[n])/y_exact[n]);
        rel_err_sym_b.push_back((y_sym_b[n]-y_exact[n])/y_exact[n]);
    }


    ofstream file;
    file.open("build/ex3_results.csv");
    file << "rel_err_norm_a,rel_err_sym_a,rel_err_norm_b,rel_err_sym_b" << endl;
    for(int n=0; n <= max_n; n++)
    {
        file << rel_err_norm_a[n] << "," << rel_err_sym_a[n] << "," << rel_err_norm_b[n] << "," << rel_err_sym_b[n] << endl;
    }
    file.close();
    
    return 0;
}