#include <iostream>
#include <eigen3/Eigen/Dense>
#include "service.cpp"

using namespace std;

int main()
{
    // a) 
    // read Image data
    Eigen::MatrixXd img;
    if(!loadData(img,"Image",512,512))
    {
       cerr << "Failed to load the File." << endl;
       return 0;
    }
    // perform SVD
    Eigen::BDCSVD<Eigen::MatrixXd> svd = img.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::VectorXd W = svd.singularValues(); // already sorted decending

    // b) 
    // rank k=(10,20,50) approximation
    vector<int> ks = {10,20,50};
    cout << "Image Compression:" << endl;
    for(int i = 0; i<ks.size(); i++)
    {
        int k = ks[i];
        cout << "k=" << k << endl;
        Eigen::MatrixXd A_k = Eigen::MatrixXd::Zero(img.cols(),img.rows());
        for(int j = 0; j<k; j++)
        {
            A_k = A_k + W(j) * U.col(j) * V.col(j).transpose();
        }
        storeData(A_k,"build/Image_"+to_string(k));
    }

    return 0;
}