#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <math.h> 


using namespace std;




double mag_mc(double H, int N)
{   
    int sig=1;// starting value of the spin is set as 1 
    int sig_tmp;
    double s_n=0;
    double E;
    double del_E = 0;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distribution(0, 1);
    std::uniform_real_distribution<double> real_dist(0,1);
    for(int i = 1 ;i<N;i++)
    {
        auto randomNumber = distribution(generator); // drawing the spin
        
        if(randomNumber == 0)
        {
            sig_tmp = -1;
        }
        else
        {
            sig_tmp = 1;
        }
        
        
        E = - sig * H;
        del_E = (-sig_tmp *H) - E;//calculating delta_E for the metropolis algorithm
        if(del_E<=0)// if the new value is energetically cheaper accept it or the new spin is the same as the old one
        {
            sig = sig_tmp;
        }
        else
        {
            auto p = real_dist(generator);//draw uniform dist. random number
            //cout<< "p  "<< p<<endl;
            if(p < exp(-del_E))// beta= 1 here , since k_B * T should = 1
            {
                sig = sig_tmp;
                
            }

        }
       
        s_n += sig; // calculating the magnetization
    }
    

   
    
    return s_n/N;
}

//function to create a linspace in c++;
vector<double> linspace(double a, double b , int N)
{   
    vector<double> linspaced;
    double delta = (b - a) / (N - 1);
    for(int i = 0 ; i < N - 1 ; i++)
    {
        linspaced.push_back(a + delta *i);

    }
    linspaced.push_back(b);

    return linspaced;
}
int main(void)
{
    
    int N = 1e5;
    vector<double> H = linspace(-5.0,5.0,1e4);//creating 1e4 values of H
    vector<double> m; // vector to store the mc_magnetization in
    auto file = fopen("build/mc_mag.dat", "w");
    fprintf(file, "# H m\n");
    for(int i = 0 ; i < 1e4 ; i++){

    m.push_back(mag_mc(H[i],N));
    fprintf(file, "%.10f %.10f\n", H[i], m[i]);
    
    }
    fclose(file);
    
    return 0;
}
