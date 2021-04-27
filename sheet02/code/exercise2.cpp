#include <iostream>
#include <vector>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "profiler.h"

using namespace std;

int main()
{
    // b) Profiling (Timing)
    int maxN = 1000;
    int numberOfTimers = 3;
    vector<vector<double>> t(numberOfTimers);
    for(int i=0; i<numberOfTimers; i++) t[i] = vector<double>(maxN);
    vector<int> Ns(maxN);
    Profiler::init(numberOfTimers);

    cout << "Starting Profiling" << endl;
    for(int N=1; N<=maxN; N++)
    {
        cout << "N=" << N << "\r";
        cout.flush();
        
        // a)
        Profiler::start(0);

        // a.1)
        Eigen::MatrixXd M = Eigen::MatrixXd::Random(N,N);
        Eigen::VectorXd b = Eigen::VectorXd::Random(N);

        Profiler::stop(0);
        Profiler::start(1);

        // a.2) LU decomposition (PartialPivLU)
        Eigen::PartialPivLU<Eigen::MatrixXd> lu = M.lu();

        Profiler::stop(1);
        Profiler::start(2);

        // a.3) solve the sle (PartialPivLU)
        Eigen::VectorXd x = lu.solve(b);

        Profiler::stop(2);

        // save the timing data
        Ns[N-1] = N;
        for(int i=0; i<numberOfTimers; i++) t[i][N] = Profiler::getTimeInS(i);
    }
    cout << endl;
    
    // write the timing data to a csv file
    ofstream file;
    file.open("build/exercise3_timing_data.csv");
    if(file.is_open())
    {
        file << "N,t_1[s],t_2[s],t_3[s]" << endl;
        for(int i=0;i<Ns.size();i++)
        {
            file << Ns[i];
            for(int j=0; j<numberOfTimers; j++)
            {
                file << "," << t[j][i];
            }
            file << endl;
        }
    }
    file.close();


    return 0;
}