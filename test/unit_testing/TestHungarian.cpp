/*******************************************************
 *
 *    _    ___  ___   ___  ___ __   __ ___  ___  ___
 *   /_\  |_ _||   \ | _ \|_ _|\ \ / /| __|| _ \/ __|
 *  / _ \  | | | |) ||   / | |  \ V / | _| |   /\__ \
 * /_/ \_\|___||___/ |_|_\|___|  \_/  |___||_|_\|___/
 *
 *
 * Copyright (C) 2025 AIOS @ AIDRIVERS Ltd - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * author = 'Adwait P. Naik'
 * email  = 'Adwait@aidrivers.ai'
 *
 *******************************************************/

#include "Hungarian.h"
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>

TEST(HungarianAlgorithmTest, Solves3x3MatrixCorrectly) {
    HungarianAlgorithm solver;

    std::vector<std::vector<double>> costMatrix = {
        {4, 1, 3},
        {2, 0, 5},
        {3, 2, 2}
    };

    std::vector<int> assignment;
    double cost = solver.Solve(costMatrix, assignment);

    // Calculate cost based on result
    double actualCost = 0;
    for (size_t i = 0; i < assignment.size(); ++i) {
        actualCost += costMatrix[i][assignment[i]];
    }

    // The expected minimal cost is 5.0, but assignment may vary due to multiple optima
    EXPECT_DOUBLE_EQ(cost, 5.0);
    EXPECT_DOUBLE_EQ(cost, actualCost);

    // Validate that the assignment is a valid permutation (i.e., each column is used once)
    std::vector<int> sortedAssignment = assignment;
    std::sort(sortedAssignment.begin(), sortedAssignment.end());
    for (size_t i = 0; i < sortedAssignment.size(); ++i) {
        EXPECT_EQ(sortedAssignment[i], i);
    }
}

// Main function for running all tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}