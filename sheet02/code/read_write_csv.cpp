#include <iostream>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <vector>

using namespace std;

// saves the given matrix to a csv file and returns whether it was successful or not
bool writeCSV(string fileName, Eigen::MatrixXd  matrix, int precision)
{
    const static Eigen::IOFormat CSVFormat(precision, Eigen::DontAlignCols, ",", "\n");
    try {
        ofstream file(fileName);
        if (file.is_open())
        {
            file << matrix.format(CSVFormat);
            file.close();
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "Writing the csv file '" << fileName << "' resulted in the Error: " << e.what() << endl;
        return false;
    }
    return true;
}

// constructs a matrix from the given csv file
Eigen::MatrixXd readCSV(string fileToOpen, bool skip_header)
{
    // inspiration:
    // https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
    // https://aleksandarhaber.com/eigen-matrix-library-c-tutorial-saving-and-loading-data-in-from-a-csv-file/

    // stores all data
    vector<double> matrixEntries;

    ifstream matrixDataFile(fileToOpen);

    string matrixRowString;
    string matrixEntry;

    int matrixRowNumber = 0;

    if(skip_header) getline(matrixDataFile, matrixRowString);
    while (getline(matrixDataFile, matrixRowString))
    {
        istringstream matrixRowStringStream(matrixRowString);
        while (getline(matrixRowStringStream, matrixEntry, ','))
        {
            matrixEntries.push_back(stod(matrixEntry));
        }
        matrixRowNumber++;
    }

    // here we convet the vector variable into the matrix and return the resulting object, 
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}