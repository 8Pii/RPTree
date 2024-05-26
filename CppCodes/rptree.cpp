#include "rptree.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

#include <map>
#include <vector>
#include <string>


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>

#include <chrono>
#include <cctype>  
#include <memory>
#include <stdexcept>
#include <functional> 


using namespace std; 
using namespace Eigen; 
using namespace utils; 
namespace fs = std::filesystem;



int globalVar_count_random = 0;
int globalVar_count_median = 0;
bool plusplus_variant = true; 


Algo_RPTree::Algo_RPTree(const std::vector<Eigen::VectorXd>& data, int minSize, double c, const std::string& angular)
    : data(data), minSize(minSize), c(c), angular(angular), err(0) {
    if (data.size() < minSize) {
        throw std::invalid_argument("Initial data size less than minSize; adjust minSize or data size.");
    }
}

// std::function<bool(const Eigen::VectorXd&)> Algo_RPTree::ChooseRule(const std::vector<Eigen::VectorXd>& S) {
//     double S_diameter = utils::diameter(S, angular);
//     double S_averageDiameter = utils::average_diameter(S, angular);
//     if (S_diameter * S_diameter <= c * S_averageDiameter * S_averageDiameter) {
//         globalVar_count_random ++;
//         Eigen::VectorXd v = Eigen::VectorXd::Random(S[0].size()).normalized(); // Random unit vector
//         std::vector<double> projections(S.size());
//         std::transform(S.begin(), S.end(), projections.begin(), [&v](const Eigen::VectorXd& vec) {
//             return vec.dot(v);
//         });
//         double median_projection = utils::median(projections);
//         return [v, median_projection](const Eigen::VectorXd& x) -> bool {
//             return x.dot(v) <= median_projection;
//         };
//     } else {
//         globalVar_count_median ++;
//         Eigen::VectorXd centroid = utils::mean(S);
//         std::vector<double> distances(S.size());
//         std::transform(S.begin(), S.end(), distances.begin(), [&centroid](const Eigen::VectorXd& vec) {
//             return (vec - centroid).norm();
//         });
//         double median_distance = utils::median(distances);
//         return [centroid, median_distance, this](const Eigen::VectorXd& x) -> bool {
//             return (x - centroid).norm() <= median_distance;
//         };
//     }
// }

std::function<bool(const Eigen::VectorXd&)> Algo_RPTree::ChooseRule(const std::vector<Eigen::VectorXd>& S) {
    double S_diameter = utils::diameter(S, angular);
    double S_averageDiameter = utils::average_diameter(S, angular);
    
    if (plusplus_variant){
        if (S_diameter * S_diameter <= c * S_averageDiameter * S_averageDiameter) {
            globalVar_count_random++;
            Eigen::VectorXd v = Eigen::VectorXd::Random(S[0].size()).normalized(); // Random unit vector
            std::vector<double> projections(S.size());
            std::transform(S.begin(), S.end(), projections.begin(), [&v](const Eigen::VectorXd& vec) {
                return vec.dot(v);
            });
            double median_projection = utils::median(projections);
            return [v, median_projection](const Eigen::VectorXd& x) -> bool {
                return x.dot(v) <= median_projection;
            };
        } else {
            globalVar_count_median++;
            Eigen::VectorXd centroid = utils::mean(S);
            std::vector<double> distances(S.size());
            std::transform(S.begin(), S.end(), distances.begin(), [&centroid](const Eigen::VectorXd& vec) {
                return (vec - centroid).norm();
            });
            double median_distance = utils::median(distances);
            return [centroid, median_distance, this](const Eigen::VectorXd& x) -> bool {
                return (x - centroid).norm() <= median_distance;
            };
        }
    } else {
        if (S_diameter * S_diameter <= c * S_averageDiameter * S_averageDiameter) {
            globalVar_count_random++;
            Eigen::VectorXd v = Eigen::VectorXd::Random(S[0].size()).normalized(); // Random unit vector
            std::vector<double> projections(S.size());
            std::transform(S.begin(), S.end(), projections.begin(), [&v](const Eigen::VectorXd& vec) {
                return vec.dot(v);
            });
            double median_projection = utils::median(projections);
            return [v, median_projection](const Eigen::VectorXd& x) -> bool {
                return x.dot(v) <= median_projection;
            };
        }  
    }

    // Throw an exception if no rule could be chosen
    throw std::runtime_error("Failed to choose a splitting rule based on the provided data.");
}

std::shared_ptr<TreeNode> Algo_RPTree::MakeTree(std::vector<Eigen::VectorXd>& S, int depth, const std::string& label) {
    if (S.empty() || S.size() < minSize) {
        return std::make_shared<TreeNode>(TreeNode{label, std::move(S), nullptr, nullptr});
    }

    auto rule = ChooseRule(S);
    std::vector<Eigen::VectorXd> S_left, S_right;

    for (auto& point : S) {
        if (rule(point)) {
            S_left.push_back(std::move(point));
        } else {
            S_right.push_back(std::move(point));
        }
    }

    auto leftChild = MakeTree(S_left, depth + 1, label + "_L");
    auto rightChild = MakeTree(S_right, depth + 1, label + "_R");

    return std::make_shared<TreeNode>(TreeNode{label, std::move(S), std::move(leftChild), std::move(rightChild)});
}

std::pair<Eigen::VectorXd, double> Algo_RPTree::compute_centroid(const std::vector<Eigen::VectorXd>& pts) {
    double QER = 0.0; 
    if (pts.empty()) {
        // std::cout << "EMPTY PT" << std::endl;
        QER = 0;
        return {Eigen::VectorXd(), QER};
    }
    
    // DEBUGGING:
    // std::cout << "PTS : ";
    // for (int i = 0; i < pts.size(); ++i) {
    //     std::cout << pts[i] << " ";
    // }
    // std::cout << std::endl;


    Eigen::VectorXd centroid = utils::mean(pts);
    for (const auto& pt : pts) {
        QER += utils::distance(centroid, pt, "NORM");
        // QER += pow((centroid - pt).squaredNorm(),0.5);   // IF ANGULAR FALSE USE THIS LINE 
        // QER += pow((centroid - pt).squaredNorm(),1);     // just to check 
        // for (int r = 0; r < pt.size(); r++){             //further checks 
        //     QER += (centroid[r] - pt[r]);
        // }
        
    }
    // DEBUGGING :
    // std::cout << "QER : " << QER << std::endl;
    // std::cout << "CENTROID VECTOR : ";
    // for (int i = 0; i < centroid.size(); ++i) {
    //     std::cout << centroid[i] << " ";
    // }
    // std::cout << std::endl;
    return {centroid, QER};
}

void Algo_RPTree::compute_error() {
    auto tree = MakeTree(this->data);   // Build the tree from the root using the entire dataset.
    process_tree(tree);                 // Recursively process each node in the tree to compute and accumulate error.
}

int Algo_RPTree::countLeafNodes(const std::shared_ptr<TreeNode>& node) {
    if (!node) return 0;
    if (node->left == nullptr && node->right == nullptr) {
        return 1;
    } else {
        return countLeafNodes(node->left) + countLeafNodes(node->right);
    }
}

void Algo_RPTree::process_tree(const std::shared_ptr<TreeNode>& node) {
    if (!node) return; // If the node is null, exit the function

    // // DEBUGGING:

    // std::cout << "NODE : ";
    // for (int i = 0; i < node.size(); ++i) {
    //     std::cout << node[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "NODE Label: " << node->label << std::endl; 

    // print elements like points (they are Eigen::VectorXd):
    // std::cout << "Points in Node: ";
    // for (const auto& point : node->points) {
    //     std::cout << point.transpose() << "; "; 
    // }
    // std::cout << std::endl;


    // Check if the node is a leaf node - A leaf node is defined as a node without any children (see linked list/array logic)
    if (node->left == nullptr && node->right == nullptr) {
        if (!node->points.empty()) {
            auto [centroid, QER] = compute_centroid(node->points);
            this->err += QER;
            this->dic[node->label] = std::make_pair(centroid, node->points); 

            // DEBUGGING 
            // std::cout << "ERR : " << this->err << std::endl;
            // std::cout << "Leaf Node Processed: " << node->label << std::endl;
            // std::cout << "Centroid: " << centroid.transpose() << std::endl;
            // std::cout << "QER: " << QER << std::endl;
            // std::cout << "Accumulated ERR: " << this->err << std::endl;
        }
    } else {
        // Recursively process the left and right children
        process_tree(node->left);
        process_tree(node->right);
    }
}

double Algo_RPTree::get_error() const {
    return this->err;
}

double cleanAndConvert(const string& input) {
    string filtered;
    for (char ch : input) {
        if (isdigit(ch) || ch == '.' || ch == '-' || ch == '\'') {
            filtered += ch;
        }
    }
    try {
        // Remove any apostrophes that could still be present
        filtered.erase(remove(filtered.begin(), filtered.end(), '\''), filtered.end());
        return stod(filtered);
    } catch (const invalid_argument& e) {
        cerr << "Conversion failed after filtering: " << filtered << " from input: " << input << endl;
        throw;  // re-throw to handle it outside if needed
    }
}




int main() {
    std::string outputPath = "/Users/OctoPii/Desktop/GMDA_RPTREE/Clustering_Results/rptree_protein_results_1000_rand.csv";
    std::ofstream resultFile(outputPath);
    resultFile << "Filename,Processing Time (min),Quantization Error,MinSize,C" << std::endl;

    // string randVsMedianPath = "/Users/OctoPii/Desktop/GMDA_RPTREE/Clustering_Results/rand_vs_median_rpt.csv";
    // ofstream randVsMedianFile(randVsMedianPath);
    // randVsMedianFile << "Filename,Random Count,Median Count" << endl;


    std::string directoryPath = "/Users/OctoPii/Desktop/GMDA_RPTREE/Protein_Data";
    // std::vector<int> minsizes = {7, 10, 15};
    // std::vector<int> c_values = {3, 10, 20};  
    std::vector<int> minsizes = {10};
    std::vector<int> c_values = {10};  

    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.path().extension() == ".csv") {
            for (int minSize : minsizes) {
                for (int c : c_values) {  // Loop through values of 'c'
                    std::string filepath = entry.path().string();
                    std::ifstream file(filepath);
                    std::vector<VectorXd> data;
                    std::string line;

                    std::cout << "Processing file: " << filepath << " with MinSize = " << minSize << " and C = " << c << std::endl;

                    if (file.is_open()) {
                        // Skip the first line (header)
                        getline(file, line);

                        while (getline(file, line)) {
                            stringstream ss(line);
                            vector<double> values;
                            string cell;

                            for (int i = 0; i < 3; ++i) {
                                if (!getline(ss, cell, ',')) {
                                    std::cerr << "Error: Incomplete line." << std::endl;
                                    return 1;
                                }
                                try {
                                    values.push_back(stod(cell));
                                } catch (const std::invalid_argument& e) {
                                    std::cerr << "Failed to convert '" << cell << "' to double: " << e.what() << std::endl;
                                    std::cerr << "Problematic line: " << line << std::endl;
                                    return 1;
                                }
                            }

                            if (!values.empty()) {
                                VectorXd vec(3);
                                for (size_t i = 0; i < values.size(); ++i) {
                                    vec(i) = values[i];
                                }
                                data.push_back(vec);
                            }
                            if (data.size() == 3000) break;
                        }
                        file.close();
                    } else {
                        std::cerr << "Failed to open file: " << filepath << std::endl;
                        continue;
                    }

                    // Display the five first lines from data
                    int count = 0;  
                    for (const auto& vec : data) {
                        if (count >= 5) break;  
                        std::cout << "Phi: " << vec(0) << ", Psi: " << vec(1) << ", Omega: " << vec(2) << std::endl;
                        count++;  
                    }

                    std::string angular = "NORM";
                    bool plusplus_variant = false; // Ensuring that its true 

                    Algo_RPTree tree(data, minSize, c, angular);
                    std::cout << "RPTree initialized. Computing error ... " << std::endl;

                    auto start = std::chrono::high_resolution_clock::now();
                    tree.compute_error();
                    auto stop = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                    double durationInMinutes = duration.count() / 60000.0;
                    double QER = tree.get_error();

                    resultFile << entry.path().filename().string() << "," 
                               << std::fixed << std::setprecision(3) << durationInMinutes << "," 
                               << QER << "," 
                               << minSize << "," 
                               << c << std::endl;

                    std::cout << "File processed: " << entry.path().filename() 
                              << " in " << durationInMinutes << " minutes with quantization error " 
                              << QER << ", MinSize " << minSize 
                              << " and C " << c << std::endl;
                    // if (plusplus_variant) {
                    //     randVsMedianFile << entry.path().filename().string() << ","
                    //                         << globalVar_count_random << ","
                    //                         << globalVar_count_median << endl;
                    //     }
                }
            }
        }
    }

    resultFile.close();
    return 0;
}