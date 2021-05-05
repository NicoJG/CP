#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Eigen>

// construction of Hamiltonian
Eigen::MatrixXd Hamiltonian(int N, double t, double eps)
{
    Eigen::MatrixXd A(N, N);

    for (int i = 0; i < N-1; i++)
    {
        A(i, i+1) = -t;
        A(i+1, i) = -t;
    }
    A(0, N-1) = -t;
    A(N-1, 0) = -t;

    A(int(N/2-1), int(N/2-1)) = eps;
    return A;
}

//Lanczos algorithm
std::tuple <Eigen::MatrixXd, Eigen::VectorXd, double> Lanczos(const Eigen::MatrixXd &A, double eps)
{   
    int k = 0;
    int N = A.cols();
    Eigen::Index groundstate;

    //T is the tridiagonal matrix, Q contains the Lanczos basis vectors q
    Eigen::MatrixXd Q  = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(N,N);

    //initialization, q_1 is a normalized random vector
    Eigen::VectorXd q_0 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd q_1 = Eigen::VectorXd::Random(N);

    q_1.normalize();
    Q.col(0) = q_1;

    double gamma = 1.;
    //define random reference-value for first calculation
    double lambda_0 = 42;

    double delta = q_1.transpose()*A*q_1;
    T(0, 0) = delta;

    for (int i = 1; i < N; i++)
    {
        Eigen::VectorXd v_1 = (A - delta*Eigen::MatrixXd::Identity(N,N))*q_1 - gamma*q_0;

        gamma = v_1.norm();

        q_0 = q_1;
        q_1 = 1./gamma*v_1;    

        Q.col(i) = q_1;

        delta = q_1.transpose()*A*q_1;

        T(i, i) = delta;
        T(i, i-1) = gamma;
        T(i-1, i) = gamma;
        
        //smallest eigenvalue and position of smallest eigenvalue
        double lambda_1 = T.block(0, 0, i, i).eigenvalues().real().minCoeff(&groundstate);

        if(std::abs(lambda_1 - lambda_0) > eps)
        {   
            lambda_0 = lambda_1;

        }
        else
        {
            k = i;
            break;
        }

    }
    
    Eigen::EigenSolver<Eigen::MatrixXd> es(T.block(0, 0, k, k));

    Eigen::VectorXd ev_groundstate = es.eigenvectors().real().col(groundstate);
    
    //couldn't find a way to calculate Q_k*ev here and return it to G.col(i) in main function
    return std::make_tuple(Q.leftCols(k), ev_groundstate, lambda_0);
}

int main()
{
    int N = 50;
    double eps = 1e-5;
    double t = 1;
    double e = -20;
    //range -20 to 20 takes 21 steps, when increment is 2
    int steps = 21;

    Eigen::MatrixXd G(N, steps);

    std::ofstream outfile;

    outfile.open("ex02_energy.log", std::fstream::out);
        outfile << "#e" << "\t" << "E_0" << "\n";
    outfile.close();
 
    for (int i = 0; i < steps; i++)
    {
        double E_0 = 0;

        Eigen::MatrixXd H = Hamiltonian(N, t, e);
    	
        Eigen::MatrixXd Q_k(N, N);
        Eigen::VectorXd ev(N);
        std::tie(Q_k, ev, E_0) = Lanczos(H, eps);
        ev.normalize();
        Eigen::VectorXd v = Q_k*ev;
        G.col(i) = v;

        outfile.open("ex02_energy.log", std::fstream::app);
            outfile << e << "\t" << E_0 << "\n";
        outfile.close();
    
        e += 2;
    }

    outfile.open("ex02_groundstate.log", std::fstream::out);
        outfile << G;
    outfile.close();

    return 0;
}