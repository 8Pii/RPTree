#include "utils.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std; 
using namespace Eigen; 



namespace utils {


double radians(double degrees) {
    return degrees * M_PI / 180.0;
}

Eigen::VectorXd mean(const std::vector<Eigen::VectorXd>& data) {
    if (data.empty()) {
        return Eigen::VectorXd();
    }
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(data[0].size());
    for (const auto& vec : data) {
        sum += vec;
    }
    Eigen::VectorXd result = sum / static_cast<double>(data.size());
    
    // Debug output
    // std::cout << "MEAN VECTOR : ";
    // for (int i = 0; i < result.size(); ++i) {
    //     std::cout << result[i] << " ";
    // }
    // std::cout << std::endl;
    
    return result;
}


double median(std::vector<double>& data) {
    double median = 0.0; 
    size_t n = data.size();
    if (n == 0) return 0; // handle empty data
    nth_element(data.begin(), data.begin() + n / 2, data.end());
    if (n % 2 == 0) {
        double mid1 = data[n / 2];
        double mid2 = *max_element(data.begin(), data.begin() + n / 2);
        median = (mid1 + mid2) / 2.0;
        // std::cout << "MEDIAN : " << median << std::endl;
        return median;
    }
    median = data[n / 2];
    // std::cout << "MEDIAN : " << median << std::endl;
    return median;
}


double distance(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const std::string& angular) {
    double sum = 0.0; 
    if (angular == "FALSE") {
        sum = (x - y).norm();
        // std::cout << "SUM (DISTANCE) : " << sum << std::endl;
        return sum;
    } else if (angular == "NORM") {
        auto sub_dist = [](double a, double b) {
            return 180.0 - std::abs(std::abs(b - a) - 180.0);
        };
        for (int i = 0; i < x.size(); ++i) {
            sum += sub_dist(x[i], y[i]);
        }
        // std::cout << "SUM (DISTANCE) : " << sum << std::endl;
        return sum;
    } else if (angular == "COS") {
        auto sub_dist = [](double a, double b) {
            return 1.0 - cos(radians(a - b));
        };
        for (int i = 0; i < x.size(); ++i) {
            sum += sub_dist(x[i], y[i]);
        }
        // std::cout << "SUM (DISTANCE) : " << sum << std::endl;
        return sum;
    } else {
        throw std::invalid_argument("Invalid angular type: must be 'FALSE', 'NORM', or 'COS'");
    }
}


double diameter(const std::vector<Eigen::VectorXd>& data, const std::string& angular) {
    double max_diam = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data.size(); ++j) {
            double dist = distance(data[i], data[j], angular);
            if (dist > max_diam) {
                max_diam = dist;
            }
        }
    }
    // std::cout << "MAX DIAMETER : " << max_diam << std::endl;
    return max_diam;
}

// double average_diameter(const vector<VectorXd>& data, const string& angular) {
//     double total_distance = 0.0;
//     size_t count = 0;

//     for (size_t i = 0; i < data.size(); ++i) {
//         for (size_t j = 0; j < data.size(); ++j) {
//             total_distance += distance(data[i], data[j], angular);
//             count++;
//         }
//     }

//     double avg_diameter = total_distance / count;
//     cout << "AVG DIAMETER : " << avg_diameter << endl;
//     return avg_diameter;
// }

double average_diameter(const vector<VectorXd>& data, const string& angular) {
    double avg_diam = 0.0;
    double dist = 0.0;
    size_t count = 0;

    for (const auto& x : data) {  // Use const auto& to avoid copying
        for (const auto& y : data) {
            dist += distance(x, y, angular);  // Ensure that utils::distance can handle Eigen::VectorXd properly
            count++;
        }
    }

    if (count > 0) {  // Protect against division by zero if data is empty
        avg_diam = dist / count;
    } else {
        avg_diam = 0.0;  // or handle this scenario as needed
    }

    // std::cout << "AVG DIAMETER : " << avg_diam << std::endl;
    return avg_diam;
}

// double average_diameter(const std::vector<Eigen::VectorXd>& data, const std::string& angular) {
//     if (data.empty()) return 0.0;
//     double total = 0.0;
//     int count = 0;
//     double avg_diam = 0.0;
//     for (size_t i = 0; i < data.size(); ++i) {
//         for (size_t j = i + 1; j < data.size(); ++j) {
//             total += distance(data[i], data[j], angular);
//             count++;
//         }
//     }
//     avg_diam = count > 0 ? total / count : 0.0;
//     std::cout << "AVG DIAM : " << avg_diam << std::endl;
//     return avg_diam;
// }


Eigen::MatrixXd compute_distance_matrix(const Eigen::MatrixXd& centroids, const std::string& angular) {
    int n = centroids.rows();
    Eigen::MatrixXd dist_matrix(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dist = distance(centroids.row(i), centroids.row(j), angular);
            dist_matrix(i, j) = dist;
            dist_matrix(j, i) = dist;
        }
    }
    return dist_matrix;
}

}


// int main() {
//     // Example usage
//     std::vector<Eigen::VectorXd> data;
//     data.reserve(10); // Reserve space for efficiency

//     // Initialize data points and push them into the vector
//     for (int i = 0; i < 10; ++i) {
//         Eigen::VectorXd point(3);
//         point << i * 3 + 1, i * 3 + 2, i * 3 + 3;
//         data.push_back(point);
//     }

//     auto avg_diam = average_diameter(data, "FALSE");
//     std::cout << "Average Diameter: " << avg_diam << std::endl;

//     return 0;
// }
