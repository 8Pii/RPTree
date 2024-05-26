#ifndef CLUSTERINGALGORITHM_HPP
#define CLUSTERINGALGORITHM_HPP

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <functional>

class ClusteringAlgorithm {
private:
    Eigen::MatrixXd data;
    int n_clusters;
    std::map<std::string, std::vector<int>> labels;
    std::map<std::string, Eigen::MatrixXd> centroids;
    std::map<std::string, std::vector<std::vector<Eigen::VectorXd>>> clusters;
    bool angular;
    std::map<std::string, std::function<void()>> algorithms;

public:
    ClusteringAlgorithm(const Eigen::MatrixXd& data, bool angular = false, int n_clusters = 0);
    void run_kmeans_plusplus();
    void run_kmeans_random();
    void run_agglomerative();
    void run_dbscan();
    void run_spectral();
    void run_gmm();
    void run_affinity();
    void organize_clusters(const std::string& algorithm_name, bool include_noise = false);
    void calculate_centroids_and_clusters();
    double compute_total_clustering_error(const std::string& algorithm_name);
    void run_algorithm(const std::string& algo_name);
};

#endif // CLUSTERINGALGORITHM_HPP
