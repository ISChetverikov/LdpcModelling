#include "../include/base_simulate.h"

std::vector<std::vector<int>> base_simulate::readAsRowSparseMatrix(std::string filename){
    std::vector<std::vector<int>> matrix;
    std::string line;
    std::ifstream myFile(filename);
    
    int val;

    while (std::getline(myFile, line))
    {
        std::stringstream ss(line);

        std::vector<int> temp;
        int columnIndex = 0;
        while (ss >> val) {
            if (val)
                temp.push_back(columnIndex);

            if (ss.peek() == ',')
                ss.ignore();

            columnIndex++;
        }

        matrix.push_back(temp);
    }

    return matrix;    
}


double base_simulate::get_sigma(double snr) {
    return sqrt(pow(10, (-snr / 10) / 2));
}


double EbN0(double snr, double R) {
    return snr - 10 * log10(R);
}