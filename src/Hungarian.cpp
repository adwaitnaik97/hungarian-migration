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

#include "Hungarian.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>

[[nodiscard]] double HungarianAlgorithm::Solve(std::vector<std::vector<double>> &DistMatrix, std::vector<int> &Assignment)
{
    std::size_t nRows = DistMatrix.size();
    std::size_t nCols = DistMatrix[0].size();
    std::vector<double> distMatrix(nRows * nCols);
    std::vector<int> assignment(nRows, -1);
    double cost = 0.0;

    for (std::size_t i = 0; i < nRows; ++i)
        for (std::size_t j = 0; j < nCols; ++j)
            distMatrix[i + nRows * j] = DistMatrix[i][j];

    assignmentOptimal(assignment, cost, distMatrix, nRows, nCols);
    Assignment = assignment;
    return cost;
}

void HungarianAlgorithm::assignmentOptimal(std::vector<int> &assignment, double &cost,
                                           const std::vector<double> &distMatrixIn,
                                           std::size_t nOfRows, std::size_t nOfColumns) const
{
    std::size_t nOfElements = nOfRows * nOfColumns;
    std::size_t minDim = std::min(nOfRows, nOfColumns);
    double minValue;

    cost = 0.0;
    assignment.assign(nOfRows, -1);
    std::vector<double> distMatrix = distMatrixIn;
    std::vector<bool> coveredColumns(nOfColumns, false), coveredRows(nOfRows, false);
    std::vector<bool> starMatrix(nOfElements, false), primeMatrix(nOfElements, false), newStarMatrix(nOfElements, false);

    if (nOfRows <= nOfColumns) {
        for (std::size_t row = 0; row < nOfRows; ++row) {
            minValue = distMatrix[row];
            for (std::size_t col = 1; col < nOfColumns; ++col)
                minValue = std::min(minValue, distMatrix[row + nOfRows * col]);
            for (std::size_t col = 0; col < nOfColumns; ++col)
                distMatrix[row + nOfRows * col] -= minValue;
        }
        for (std::size_t row = 0; row < nOfRows; ++row)
            for (std::size_t col = 0; col < nOfColumns; ++col)
                if (std::fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON && !coveredColumns[col]) {
                    starMatrix[row + nOfRows * col] = true;
                    coveredColumns[col] = true;
                    break;
                }
    } else {
        for (std::size_t col = 0; col < nOfColumns; ++col) {
            minValue = distMatrix[nOfRows * col];
            for (std::size_t row = 1; row < nOfRows; ++row)
                minValue = std::min(minValue, distMatrix[row + nOfRows * col]);
            for (std::size_t row = 0; row < nOfRows; ++row)
                distMatrix[row + nOfRows * col] -= minValue;
        }
        for (std::size_t col = 0; col < nOfColumns; ++col)
            for (std::size_t row = 0; row < nOfRows; ++row)
                if (std::fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON && !coveredRows[row]) {
                    starMatrix[row + nOfRows * col] = true;
                    coveredRows[row] = true;
                    coveredColumns[col] = true;
                    break;
                }
        std::fill(coveredRows.begin(), coveredRows.end(), false);
    }

    step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
           coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

    computeAssignmentCost(assignment, cost, distMatrixIn, nOfRows);
}

void HungarianAlgorithm::buildAssignmentVector(std::vector<int> &assignment, const std::vector<bool> &starMatrix,
                                               std::size_t nOfRows, std::size_t nOfColumns) const
{
    for (std::size_t row = 0; row < nOfRows; ++row)
        for (std::size_t col = 0; col < nOfColumns; ++col)
            if (starMatrix[row + nOfRows * col]) {
                assignment[row] = static_cast<int>(col);
                break;
            }
}

void HungarianAlgorithm::computeAssignmentCost(const std::vector<int> &assignment, double &cost,
                                               const std::vector<double> &distMatrix, std::size_t nOfRows) const
{
    cost = 0.0;
    for (std::size_t row = 0; row < nOfRows; ++row) {
        int col = assignment[row];
        if (col >= 0)
            cost += distMatrix[row + nOfRows * col];
    }
}

void HungarianAlgorithm::step2a(std::vector<int> &assignment, std::vector<double> &distMatrix,
                                std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                                std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                                std::vector<bool> &coveredRows, std::size_t nOfRows,
                                std::size_t nOfColumns, std::size_t minDim) const
{
    for (std::size_t col = 0; col < nOfColumns; ++col)
        for (std::size_t row = 0; row < nOfRows; ++row)
            if (starMatrix[row + nOfRows * col]) {
                coveredColumns[col] = true;
                break;
            }

    step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
           coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

void HungarianAlgorithm::step2b(std::vector<int> &assignment, std::vector<double> &distMatrix,
                                std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                                std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                                std::vector<bool> &coveredRows, std::size_t nOfRows,
                                std::size_t nOfColumns, std::size_t minDim) const
{
    std::size_t coveredCount = std::count(coveredColumns.begin(), coveredColumns.end(), true);
    if (coveredCount == minDim) {
        buildAssignmentVector(assignment, starMatrix, nOfRows, nOfColumns);
        return;
    }
    step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
          coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

void HungarianAlgorithm::step3(std::vector<int> &assignment, std::vector<double> &distMatrix,
                               std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                               std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                               std::vector<bool> &coveredRows, std::size_t nOfRows,
                               std::size_t nOfColumns, std::size_t minDim) const
{
    while (true) {
        bool found = false;
        for (std::size_t col = 0; col < nOfColumns; ++col) {
            if (!coveredColumns[col]) {
                for (std::size_t row = 0; row < nOfRows; ++row) {
                    if (!coveredRows[row] && std::fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON) {
                        primeMatrix[row + nOfRows * col] = true;
                        std::size_t starCol;
                        for (starCol = 0; starCol < nOfColumns; ++starCol)
                            if (starMatrix[row + nOfRows * starCol]) break;
                        if (starCol == nOfColumns) {
                            step4(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
                                  coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
                            return;
                        } else {
                            coveredRows[row] = true;
                            coveredColumns[starCol] = false;
                            found = true;
                            break;
                        }
                    }
                }
            }
        }
        if (!found) break;
    }
    step5(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
          coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

void HungarianAlgorithm::step4(std::vector<int> &assignment, std::vector<double> &distMatrix,
                               std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                               std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                               std::vector<bool> &coveredRows, std::size_t nOfRows,
                               std::size_t nOfColumns, std::size_t minDim,
                               std::size_t row, std::size_t col) const
{
    newStarMatrix = starMatrix;
    newStarMatrix[row + nOfRows * col] = true;
    std::size_t starCol = col, starRow = 0;

    while (true) {
        for (starRow = 0; starRow < nOfRows; ++starRow)
            if (starMatrix[starRow + nOfRows * starCol]) break;
        if (starRow == nOfRows) break;
        newStarMatrix[starRow + nOfRows * starCol] = false;
        std::size_t primeCol = 0;
        for (; primeCol < nOfColumns; ++primeCol)
            if (primeMatrix[starRow + nOfRows * primeCol]) break;
        newStarMatrix[starRow + nOfRows * primeCol] = true;
        starCol = primeCol;
    }

    starMatrix = newStarMatrix;
    std::fill(primeMatrix.begin(), primeMatrix.end(), false);
    std::fill(coveredRows.begin(), coveredRows.end(), false);
    std::fill(coveredColumns.begin(), coveredColumns.end(), false);

    step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
           coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

void HungarianAlgorithm::step5(std::vector<int> &assignment, std::vector<double> &distMatrix,
                               std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                               std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                               std::vector<bool> &coveredRows, std::size_t nOfRows,
                               std::size_t nOfColumns, std::size_t minDim) const
{
    double minUncovered = DBL_MAX;
    for (std::size_t row = 0; row < nOfRows; ++row)
        if (!coveredRows[row])
            for (std::size_t col = 0; col < nOfColumns; ++col)
                if (!coveredColumns[col])
                    minUncovered = std::min(minUncovered, distMatrix[row + nOfRows * col]);

    for (std::size_t row = 0; row < nOfRows; ++row)
        for (std::size_t col = 0; col < nOfColumns; ++col) {
            if (coveredRows[row]) distMatrix[row + nOfRows * col] += minUncovered;
            if (!coveredColumns[col]) distMatrix[row + nOfRows * col] -= minUncovered;
        }

    step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
          coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}