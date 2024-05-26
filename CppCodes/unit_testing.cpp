#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "utils.hpp"  // Assuming this is the location of your utility functions

using namespace std;
using namespace Eigen;
using namespace utils;

// Function to compare two doubles with a tolerance
bool compareDouble(double x, double y, double tol = 1e-6) {
    return fabs(x - y) < tol;
}

// Function to run all the tests
void runTests() {
    // Test data setup
    vector<VectorXd> data;
    VectorXd v1(3), v2(3), v3(3);
    v1 << 1, 2, 3;
    v2 << 4, 5, 6;
    v3 << 7, 8, 9;
    data.push_back(v1);
    data.push_back(v2);
    data.push_back(v3);

    // Test Mean
    VectorXd meanResult = mean(data);
    VectorXd expectedMean(3);
    expectedMean << 4, 5, 6;
    cout << "Test Mean: " << (meanResult.isApprox(expectedMean) ? "PASS" : "FAIL") << endl;

    // Test Median
    vector<double> medData = {1, 2, 3, 4, 5};
    double medianResult = median(medData);
    double expectedMedian = 3;
    cout << "Test Median: " << (compareDouble(medianResult, expectedMedian) ? "PASS" : "FAIL") << endl;

    // Test Distance
    double distanceResult = distance(v1, v2, "FALSE");
    double expectedDistance = sqrt(27);
    cout << "Test Distance: " << (compareDouble(distanceResult, expectedDistance) ? "PASS" : "FAIL") << endl;

    // Test Diameter
    double diameterResult = diameter(data, "FALSE");
    double expectedDiameter = sqrt(108);
    cout << "Test Diameter: " << (compareDouble(diameterResult, expectedDiameter) ? "PASS" : "FAIL") << endl;

    // Test Average Diameter
    double avgDiameterResult = average_diameter(data, "FALSE");
    double expectedAvgDiameter = 4.618802153517006;//(sqrt(27) + sqrt(108) + sqrt(27)) / 3;
    cout << "Test Average Diameter: " << (compareDouble(avgDiameterResult, expectedAvgDiameter) ? "PASS" : "FAIL") << endl;
}

int main() {
    runTests();  // Execute all the tests
    return 0;
}
