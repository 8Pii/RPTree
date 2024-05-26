//utils.hpp

#ifndef UTILS_HPP
#define UTILS_HPP


#include <Eigen/Dense>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>


namespace utils{

double radians(double degrees);    

Eigen::VectorXd mean(const std::vector<Eigen::VectorXd>& data);

double median(std::vector<double>& data);

double distance(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const std::string& angular);

double diameter(const std::vector<Eigen::VectorXd>& data, const std::string& angular);

// double average_diameter(const std::vector<Eigen::VectorXd>& data, const std::string& angular);

double average_diameter(const std::vector<Eigen::VectorXd>& data, const std::string& angular);

Eigen::MatrixXd compute_distance_matrix(const Eigen::MatrixXd& centroids, const std::string& angular_type);
}


#endif // UTILS_HPP