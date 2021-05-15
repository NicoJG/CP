#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <fstream>

const std::complex<double> i_c(0., 1.);

//reverse binary number
int reverseNum(int num, int m)
{
    int a[m];
    for (int i = 0; i < m; i++)
    {
        a[i] = num  % 2;
        num = num/2;
    }

    int num_reverse = 0;
    int i_reverse = m-1;
    for(int i = 0; i < m; i++)
    {   
        num_reverse += pow(2, i)*a[i_reverse];
        i_reverse--;
    }
    return num_reverse;
}

//function f
double f(int l)
{
    return sqrt(1 + l);
}

//discrete FFT
Eigen::VectorXcd dFFT(double (* func)(int), int m)
{       
    //define S_0
    int N = pow(2, m);
    Eigen::MatrixXcd S_0(N, N);
 
    for(int l = 0; l < N; l++)
    {
        double tmp = f(reverseNum(l, m));
        for(int j = 0; j < N; j++)
            {
                S_0(j, l) = tmp;
            }
    }
    //define S_k
    for (int k = 1; k <= m; k++)
    {
        Eigen::MatrixXcd S_k(N, int(pow(2, m-k)));
        for (int j = 0; j < pow(2,k); j++)
        {
            for (int l = 0; l < pow(2, m-k); l++)
            {
                S_k(j, l) = S_0(j, 2*l) + exp(2 * M_PI * i_c * double(j) / pow(2, k)) * S_0(j, 2*l+1);
                for(int i = 1; i < pow(2, m-k); i++)
                //inner loop for j >= 2^k - 1
                {
                    S_k(j + i*pow(2,k), l) = S_k(j, l);
                }
            }
        }
        S_0 = S_k;
    }
    //k = m
    Eigen::VectorXcd F = S_0.col(0);
    return F;
}

//direct evaluation of discrete FT
Eigen::VectorXcd dFT(double (* func)(int), int m)
{
    int N = pow(2, m);
    Eigen::VectorXcd F = Eigen::VectorXcd::Zero(N);
    for (int j = 0; j < N; j++)
    {
        for (int l = 0; l < N; l++)
        {
            F(j) +=  exp(2*M_PI*i_c*double(j)*double(l)/double(N))*func(l);
        }
    }
    return F;
}

int main()
{    
    std::ofstream outfile;

    for (int m = 3; m <= 4; m++)
    {
        Eigen::VectorXcd F_fft = dFFT(f, m);
        Eigen::VectorXcd F_ft = dFT(f, m);

        outfile.open("build/ex02_m" + std::to_string(m) + ".log", std::fstream::out);
        outfile << "#R(FFT)"<< "\t" << "Im(FFT)" << "\t" << "R(FT)" << "\t" << "Im(FT)" <<"\n" ;            

        for (int i = 0; i < F_fft.size(); i++)
        {
            outfile << F_fft(i).real() << "\t" << F_fft(i).imag()  << "\t" << F_ft(i).real() << "\t" << F_ft(i).imag() <<"\n" ;            
        }
        outfile.close();
    }
    return 0;
}