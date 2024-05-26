#include "kmeanspp.hpp"

//Data structures & functions
#include <map>
#include <cmath>
#include <vector>
#include <cctype>  
#include <string>
#include <cstdlib> 
#include <algorithm>
#include <functional> 

//Input / Output manip
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>

//System
#include <chrono>
#include <memory>
#include <stdexcept>
#include <filesystem>



using namespace std; 
namespace fs = std::filesystem;


KmeansPlusPlus::KmeansPlusPlus(int k, int pointnumber) : Kmeans(k, pointnumber) {}

KmeansPlusPlus::KmeansPlusPlus(int k) : Kmeans(k) {}

int KmeansPlusPlus::NearestCenter(Point &p, int alreadyInitCenterNumber, float &minDistance) {
    minDistance = std::numeric_limits<float>::max();
    int k_id = -1;
    float dis;
    auto centersIter = _Centers.begin();
    std::advance(centersIter, alreadyInitCenterNumber);  // Adjust to start from the correct position
    for (int k = 0; k <= alreadyInitCenterNumber && centersIter != _Centers.end(); ++centersIter, ++k) {
        dis = Distance(p, *centersIter);
        if (dis < minDistance) {
            minDistance = dis;
            k_id = k;
        }
    }
    return k_id;
}

void KmeansPlusPlus::InitCenters() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(0, _PointNumber - 1);

    // Select the first center randomly
    int id = distribution(gen);
    auto it = _Points.begin();
    std::advance(it, id);  // Move iterator to the randomly selected starting point
    _Centers.push_back(*it);

    // Prepare for choosing subsequent centers
    for (int i = 1; i < _K; ++i) {
        float sum = 0;
        std::vector<float> distances(_PointNumber);
        auto pointIter = _Points.begin();

        // Compute distances from each point to nearest center
        for (int j = 0; j < _PointNumber && pointIter != _Points.end(); ++pointIter, ++j) {
            float min_distance;
            NearestCenter(*pointIter, i - 1, min_distance);
            distances[j] = min_distance * min_distance;  // square the distances to bias towards larger gaps
            sum += distances[j];
        }

        // Select a new center based on weighted probabilities
        std::uniform_real_distribution<float> dist(0, sum);
        float target = dist(gen);
        float cumulative = 0;
        pointIter = _Points.begin();
        for (int j = 0; j < _PointNumber && pointIter != _Points.end(); ++pointIter, ++j) {
            cumulative += distances[j];
            if (cumulative >= target) {
                _Centers.push_back(*pointIter);
                break;
            }
        }
    }
}



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " output_base_name" << std::endl;
        return 1;
    }

    std::string outputBaseName = argv[1];
    std::string outputPath = "/Users/OctoPii/Desktop/GMDA_RPTREE/Clustering_Results/" + outputBaseName + "_protein_results_10k.csv";

    std::ofstream resultFile(outputPath);
    resultFile << "Filename,Processing Time (min),Clustering Error,MINSIZE" << std::endl;

    std::string directoryPath = "/Users/OctoPii/Desktop/GMDA_RPTREE/Protein_Data";

    std::vector<int> minsizes = {7, 10, 15};


    std::cout << "Output file : " << outputBaseName << "_protein_results.csv created successfully." << std::endl;
    std::cout << "STARTING READ AND PROCESSING DATA." << std::endl;
    std::cout << "\n" << std::endl;

    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.path().extension() == ".csv") {
            std::string filepath = entry.path().string();

            for (int MINSIZE : minsizes) {
                std::ifstream file(filepath);
                std::list<Point> data;
                std::string line;

                std::cout << "Processing file: " << filepath << " with MINSIZE = " << MINSIZE << std::endl;

                if (file.is_open()) {
                    getline(file, line);  // Skip the first line (header)

                    while (getline(file, line)) {
                        std::stringstream ss(line);
                        std::vector<double> values;
                        std::string cell;

                        for (int i = 0; i < 3; ++i) {
                            if (!getline(ss, cell, ',')) {
                                std::cerr << "Error: Incomplete line." << std::endl;
                                return 1;
                            }

                            try {
                                values.push_back(std::stod(cell));
                            } catch (const std::invalid_argument& e) {
                                std::cerr << "Failed to convert '" << cell << "' to double: " << e.what() << std::endl;
                                std::cerr << "Problematic line: " << line << std::endl;
                                return 1;
                            }
                        }

                        if (!values.empty()) {
                            Point p(values[0], values[1], values[2]);
                            data.push_back(p);
                        } 
                        if (data.size() == 10000) break;        //for comparision with rptree
                    } 
                    file.close();
                } else {
                    std::cerr << "Failed to open file: " << filepath << std::endl;
                    continue;
                }

                // Display the first five points from data
                auto it = data.begin();
                int count = 0;
                for (; it != data.end() && count < 5; ++it, ++count) {
                    cout << "Phi: " << (*it)[0] << ", Psi: " << (*it)[1] << ", Omega: " << (*it)[2] << endl;
                }
                
                int DEPTH_MAX = int((log(data.size() / MINSIZE)) / log(2)) + 1;
                int N_CLUSTER_MAX = pow(DEPTH_MAX, 2);
                int k = N_CLUSTER_MAX;

                Kmeans* kmeansPtr;

                if (outputBaseName == "plusplus") {
                    KmeansPlusPlus kmeans(k, data.size());
                    kmeansPtr = &kmeans;
                } else {
                    Kmeans kmeans(k, data.size());
                    kmeansPtr = &kmeans;
                }

                kmeansPtr->InitPoints(data); // Initialize Kmeans with the data
                kmeansPtr->InitCenters(); // Initialize centers
                std::cout << "Kmeans initialized. Computing error ..." << std::endl;

                auto start = std::chrono::high_resolution_clock::now();
                kmeansPtr->RunKmean();
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                double durationInMinutes = duration.count() / 60000.0;

                double totalError = kmeansPtr->CalculateTotalClusteringError();
                resultFile << entry.path().filename().string() << "," << std::fixed << std::setprecision(3) << durationInMinutes << "," << totalError << "," << MINSIZE << std::endl;

                std::cout << "File processed: " << entry.path().filename() << " in " << durationInMinutes << " minutes with error " << totalError << " and MINSIZE " << MINSIZE << std::endl;
                std::cout << "\n" << std::endl;
            }
        }
    }
    resultFile.close();
    return 0;
}