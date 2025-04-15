///////////////////////////////////////////////////////////////////////////////
// Hungarian.cpp: Implementation file for Class HungarianAlgorithm.
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

#include <cstdlib> // For standard library functions
#include <cfloat>  // For DBL_MAX
#include <cmath>   // For mathematical functions

#include <algorithm>

//********************************************************//
// A single function wrapper for solving assignment problem.
//********************************************************//
double HungarianAlgorithm::Solve(std::vector<std::vector<double>> &DistMatrix, std::vector<int> &Assignment)
{
    double cost = 0.0;

    std::size_t nRows = DistMatrix.size();
    std::size_t nCols = DistMatrix[0].size();

    std::vector<double> distMatrixIn(nRows * nCols);
    std::vector<int> assignment(nRows);

    // Fill the distMatrixIn. Mind the index is "i + j * nRows".
    for (std::size_t i = 0; i < nRows; i++)
    {
        for (std::size_t j = 0; j < nCols; j++)
        {
            distMatrixIn[i + j * nRows] = DistMatrix[i][j];
        }
    }

    // Call solving function
    // Since the raw pointers (int *assignment, double *distMatrix) are rewritten as std::vector
    // The memory management is done automatically. Hence, we don't need new and delete.
    assignmentOptimal(assignment, cost, distMatrixIn, nRows, nCols);

    Assignment.clear();
    for (std::size_t r = 0; r < nRows; r++)
    {
        Assignment.push_back(assignment[r]);
    }

    return cost;
}

//********************************************************//
// Solve optimal solution for assignment problem using Munkres algorithm, also
// known as Hungarian Algorithm.
//********************************************************//
void HungarianAlgorithm::assignmentOptimal(std::vector<int> &assignment, double &cost,
                                           const std::vector<double> &distMatrixIn, std::size_t nOfRows, std::size_t nOfColumns) const
{
    std::size_t nOfElements = nOfRows * nOfColumns;
    std::size_t minDim;
    double value, minValue;

    /* Initialization */
    cost = 0.0;
    assignment.assign(nOfRows, -1);

    /* Generate working copy of distance matrix */
    std::vector<double> distMatrix(nOfElements);
    for (std::size_t i = 0; i < nOfElements; i++)
    {
        value = distMatrixIn[i];
        if (value < 0)
        {
            std::cerr << "All matrix elements have to be non-negative." << std::endl;
        }
        distMatrix[i] = value;
    }

    /* Memory allocation */
    std::vector<bool> coveredColumns(nOfColumns, false);
    std::vector<bool> coveredRows(nOfRows, false);
    std::vector<bool> starMatrix(nOfElements, false);
    std::vector<bool> primeMatrix(nOfElements, false);
    std::vector<bool> newStarMatrix(nOfElements, false); // Used in step 4

    /* Preliminary steps */
    if (nOfRows <= nOfColumns)
    {
        minDim = nOfRows;

        for (std::size_t row = 0; row < nOfRows; row++)
        {
            // Find the smallest element in the row
            minValue = distMatrix[row];
            for (std::size_t col = 1; col < nOfColumns; col++)
            {
                value = distMatrix[row + nOfRows * col];
                if (value < minValue)
                {
                    minValue = value;
                }
            }

            // Subtract the smallest element from each element of the row
            for (std::size_t col = 0; col < nOfColumns; col++)
            {
                distMatrix[row + nOfRows * col] -= minValue;
            }
        }

        // Steps 1 and 2a
        for (std::size_t row = 0; row < nOfRows; row++)
        {
            for (std::size_t col = 0; col < nOfColumns; col++)
            {
                if (std::fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON && !coveredColumns[col])
                {
                    starMatrix[row + nOfRows * col] = true;
                    coveredColumns[col] = true;
                    break;
                }
            }
        }
    }
    else
    {
        minDim = nOfColumns;

        for (std::size_t col = 0; col < nOfColumns; col++)
        {
            // Find the smallest element in the column
            minValue = distMatrix[nOfRows * col];
            for (std::size_t row = 1; row < nOfRows; row++)
            {
                value = distMatrix[row + nOfRows * col];
                if (value < minValue)
                {
                    minValue = value;
                }
            }

            // Subtract the smallest element from each element of the column
            for (std::size_t row = 0; row < nOfRows; row++)
            {
                distMatrix[row + nOfRows * col] -= minValue;
            }
        }

        // Steps 1 and 2a
        for (std::size_t col = 0; col < nOfColumns; col++)
        {
            for (std::size_t row = 0; row < nOfRows; row++)
            {
                if (std::fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON && !coveredRows[row])
                {
                    starMatrix[row + nOfRows * col] = true;
                    coveredColumns[col] = true;
                    coveredRows[row] = true;
                    break;
                }
            }
        }

        std::fill(coveredRows.begin(), coveredRows.end(), false);
    }

    /* Move to step 2b */
    step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
           coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);

    /* Compute cost and remove invalid assignments */
    computeAssignmentCost(assignment, cost, distMatrixIn, nOfRows);
}

/********************************************************/
void HungarianAlgorithm::buildAssignmentVector(std::vector<int> &assignment,
                                               const std::vector<bool> &starMatrix,
                                               std::size_t nOfRows,
                                               std::size_t nOfColumns) const
{
    for (std::size_t row = 0; row < nOfRows; row++)
    {
        for (std::size_t col = 0; col < nOfColumns; col++)
        {
            if (starMatrix[row + nOfRows * col])
            {
#ifdef ONE_INDEXING
                assignment[row] = col + 1; // MATLAB-Indexing
#else
                assignment[row] = col;
#endif
                break;
            }
        }
    }
}

/********************************************************/
void HungarianAlgorithm::computeAssignmentCost(const std::vector<int> &assignment, double &cost,
                                               const std::vector<double> &distMatrix,
                                               std::size_t nOfRows) const
{
    cost = 0.0;
    for (std::size_t row = 0; row < nOfRows; row++)
    {
        std::size_t col = assignment[row];
        if (col >= 0)
        {
            cost += distMatrix[row + nOfRows * col];
        }
    }
}

/********************************************************/
void HungarianAlgorithm::step2a(std::vector<int> &assignment,
                                std::vector<double> &distMatrix,
                                std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                                std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                                std::vector<bool> &coveredRows, std::size_t nOfRows,
                                std::size_t nOfColumns, std::size_t minDim) const
{

    // std::vector<bool> starMatrixTemp;
    // std::vector<bool> columnEnd;

    /* cover every column containing a starred zero */
    for (std::size_t col = 0; col < nOfColumns; col++)
    {
        for (std::size_t row = 0; row < nOfRows; row++)
        {
            if (starMatrix[row + nOfRows * col])
            {
                coveredColumns[col] = true;
                break;
            }
        }
    }

    /* move to step 3 */
    step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
           coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void HungarianAlgorithm::step2b(std::vector<int> &assignment,
                                std::vector<double> &distMatrix,
                                std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                                std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                                std::vector<bool> &coveredRows, std::size_t nOfRows,
                                std::size_t nOfColumns, std::size_t minDim) const
{

    /* count covered columns */
    std::size_t nOfCoveredColumns = 0;

    for (std::size_t col = 0; col < nOfColumns; col++)
    {
        if (coveredColumns[col])
        {
            nOfCoveredColumns++;
        }
    }

    if (nOfCoveredColumns == minDim)
    {
        /* algorithm finished */
        buildAssignmentVector(assignment, starMatrix, nOfRows, nOfColumns);
        return;
    }
    else
    {
        /* move to step 3 */
        step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
              coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
    }
}

/********************************************************/
void HungarianAlgorithm::step3(std::vector<int> &assignment,
                               std::vector<double> &distMatrix,
                               std::vector<bool> &starMatrix, std::vector<bool> &newStarMatrix,
                               std::vector<bool> &primeMatrix, std::vector<bool> &coveredColumns,
                               std::vector<bool> &coveredRows, std::size_t nOfRows,
                               std::size_t nOfColumns, std::size_t minDim) const
{

    bool zerosFound = true;
    std::size_t row, col, starCol;

    while (zerosFound)
    {
        zerosFound = false;
        for (std::size_t col = 0; col < nOfColumns; col++)
        {
            if (!coveredColumns[col])
            {
                for (std::size_t row = 0; row < nOfRows; row++)
                {
                    if ((!coveredRows[row]) && (fabs(distMatrix[row + nOfRows * col]) < DBL_EPSILON))
                    {
                        primeMatrix[row + nOfRows * col] = true;

                        /* find starred zero in current row */
                        for (std::size_t starCol = 0; starCol < nOfColumns; starCol++)
                        {
                            if (starMatrix[row + nOfRows * starCol])
                                break;
                        }

                        if (starCol == nOfColumns) /* no starred zero found */
                        {
                            /* move to step 4 */
                            step4(assignment, distMatrix, starMatrix, newStarMatrix,
                                  primeMatrix, coveredColumns, coveredRows, nOfRows,
                                  nOfColumns, minDim, row, col);
                            return;
                        }
                        else
                        {
                            /* cover row and uncover the star of the column */
                            coveredRows[row] = true;
                            coveredColumns[starCol] = false;
                            zerosFound = true;
                            break; // Break out of the inner loop
                        }
                    }
                }
            }
        }
    }
}

/********************************************************/
void HungarianAlgorithm::step4(std::vector<int> &assignment,
                               std::vector<double> &distMatrix,
                               std::vector<bool> &starMatrix,
                               std::vector<bool> &newStarMatrix,
                               std::vector<bool> &primeMatrix,
                               std::vector<bool> &coveredColumns,
                               std::vector<bool> &coveredRows,
                               std::size_t nOfRows, std::size_t nOfColumns,
                               std::size_t minDim, std::size_t row,
                               std::size_t col) const
{
    std::size_t nOfElements = nOfRows * nOfColumns;

    /* generate temporary copy of starMatrix */
    newStarMatrix = starMatrix;

    /* star current zero */
    newStarMatrix[row + nOfRows * col] = true;

    /* find starred zero in current column */
    std::size_t starCol = col;
    std::size_t starRow;
    for (starRow = 0; starRow < nOfRows; ++starRow)
    {
        if (starMatrix[starRow + nOfRows * starCol])
        {
            break;
        }
    }

    while (starRow < nOfRows)
    {
        /* unstar the starred zero */
        newStarMatrix[starRow + nOfRows * starCol] = false;

        /* find primed zero in current row */
        std::size_t primeRow = starRow;
        std::size_t primeCol;
        for (primeCol = 0; primeCol < nOfColumns; ++primeCol)
        {
            if (primeMatrix[primeRow + nOfRows * primeCol])
            {
                break;
            }
        }

        /* star the primed zero */
        newStarMatrix[primeRow + nOfRows * primeCol] = true;

        /* find starred zero in current column */
        starCol = primeCol;
        for (starRow = 0; starRow < nOfRows; ++starRow)
        {
            if (starMatrix[starRow + nOfRows * starCol])
            {
                break;
            }
        }
    }

    /* use temporary copy as new starMatrix */
    /* delete all primes, uncover all rows */
    std::fill(primeMatrix.begin(), primeMatrix.end(), false);
    starMatrix = newStarMatrix;
    std::fill(coveredRows.begin(), coveredRows.end(), false);

    /* move to step 2a */
    step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
           coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

/********************************************************/
void HungarianAlgorithm::step5(std::vector<int> &assignment,
                               std::vector<double> &distMatrix,
                               std::vector<bool> &starMatrix,
                               std::vector<bool> &newStarMatrix,
                               std::vector<bool> &primeMatrix,
                               std::vector<bool> &coveredColumns,
                               std::vector<bool> &coveredRows,
                               std::size_t nOfRows, std::size_t nOfColumns,
                               std::size_t minDim) const
{
    double h = DBL_MAX;
    for (std::size_t row = 0; row < nOfRows; row++)
    {
        if (!coveredRows[row])
        {
            for (std::size_t col = 0; col < nOfColumns; col++)
            {
                if (!coveredColumns[col])
                {
                    h = std::min(h, distMatrix[row + nOfRows * col]);
                }
            }
        }
    }

    /* add h to every element of covered rows */
    for (std::size_t row = 0; row < nOfRows; row++)
    {
        if (coveredRows[row])
        {
            for (std::size_t col = 0; col < nOfColumns; col++)
            {
                distMatrix[row + nOfRows * col] += h;
            }
        }
    }

    /* subtract h from every element of uncovered columns */
    for (std::size_t col = 0; col < nOfColumns; col++)
    {
        if (!coveredColumns[col])
        {
            for (std::size_t row = 0; row < nOfRows; row++)
            {
                distMatrix[row + nOfRows * col] -= h;
            }
        }
    }

    /* repeat step 3 */
    step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix,
          coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}