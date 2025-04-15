///////////////////////////////////////////////////////////////////////////////
// Hungarian.h: Header file for Class HungarianAlgorithm.
//
// This is a C++ wrapper with slight modification of a hungarian algorithm
// implementation by Markus Buehren. The original implementation is a few
// mex-functions for use in MATLAB, found here:
// http://www.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem
//
// Both this code and the orignal code are published under the BSD license.
// by Cong Ma, 2016
//
// modified by Adwait Naik, AI Drivers 2025
//

#ifndef __HUNGARIAN_ALGORITHM_H_
#define __HUNGARIAN_ALGORITHM_H_

#include <vector>
#include <iostream>
#include <limits>
#include <memory>

class HungarianAlgorithm
{

public:
    HungarianAlgorithm() = default;

    ~HungarianAlgorithm() = default;

    // Solve the assignment problem
    double Solve(std::vector<std::vector<double>> &DistMatrix, std::vector<int> &Assignment);

private:
    void assignmentOptimal(std::vector<int> &assignment, double &cost, const std::vector<double> &distMatrix, std::size_t nRows, std::size_t nCols) const;

    void buildAssignmentVector(std::vector<int> &assignment, const std::vector<bool> &starMatrix, std::size_t nOfRows, std::size_t nOfCols) const;

    void computeAssignmentCost(const std::vector<int> &assignment, double &cost,
                               const std::vector<double> &distMatrix,
                               std::size_t nOfRows) const;

    void step2a(std::vector<int> &assignment,
                std::vector<double> &distMatrix,
                std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                std::vector<bool> &coveredRows, std::size_t nOfRows,
                std::size_t nOfColumns, std::size_t minDim) const;

    void step2b(std::vector<int> &assignment,
                std::vector<double> &distMatrix,
                std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                std::vector<bool> &coveredRows, std::size_t nOfRows,
                std::size_t nOfColumns, std::size_t minDim) const;

    void step3(std::vector<int> &assignment,
               std::vector<double> &distMatrix,
               std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
               std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
               std::vector<bool> &coveredRows, std::size_t nOfRows,
               std::size_t nOfColumns, std::size_t minDim) const;

    void step4(std::vector<int> &assignment,
               std::vector<double> &distMatrix,
               std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
               std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
               std::vector<bool> &coveredRows, std::size_t nOfRows,
               std::size_t nOfColumns, std::size_t minDim, std::size_t row,
               std::size_t col) const;

    void step5(std::vector<int> &assignment,
               std::vector<double> &distMatrix,
               std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
               std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
               std::vector<bool> &coveredRows, std::size_t nOfRows,
               std::size_t nOfColumns, std::size_t minDim) const;
};

#endif // __HUNGARIAN_ALGORITHM_H_
