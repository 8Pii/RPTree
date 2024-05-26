#include "kmeans.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cctype>  
#include <memory>
#include <stdexcept>
#include <functional> 

using namespace std; 


Kmeans::Kmeans(int k, int pointnumber) : _K(k), _PointNumber(pointnumber) {}

struct ClusteringMetrics {
    float totalWithinClusterError;
    float totalBetweenClusterDistance;
};

float Kmeans::Distance(const Point &p1, const Point &p2)
{
    auto sub_dist = [](float a, float b){
        return 180 - std::abs(std::abs(b-a)-180);
    };
    float sum = 0.0; 
    for (int i = 0; i < 3; i++) {  
        sum += sub_dist(p1[i], p2[i]);
    }
    return sum;
}

void Kmeans::InitPoints() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0, 100.0);
    for (int i = 0; i < _PointNumber; i++) {
        Point p(distribution(gen), distribution(gen), distribution(gen));
        _Points.push_back(p);
    }
}

void Kmeans::InitPoints(const std::list<Point>& pointlist) {
    _Points = pointlist;
    _PointNumber = pointlist.size();
}

void Kmeans::InitPoints(std::string fileName) {
    std::ifstream fs(fileName.c_str(), std::ifstream::in);
    if (fs.is_open()) {
        float x, y, z;
        while (fs >> x >> y >> z) {
            _Points.push_back(Point(x, y, z));
        }
        _PointNumber = _Points.size();
        fs.close();
    } else {
        std::cout << fileName << " error opening file\n";
    }
}

void Kmeans::InitCenters() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(0, _PointNumber - 1);
    std::map<int, int> uniqueMap;
    while (uniqueMap.size() < _K) {
        int id = distribution(gen);
        uniqueMap[id] = id;
    }
    for (auto id : uniqueMap) {
        auto it = _Points.begin();
        std::advance(it, id.first);
        _Centers.push_back(*it);
    }
}

void Kmeans::InitSpecifiedCenters() {
#ifdef USING_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < _K; i++) {
        auto it = _Points.begin();
        std::advance(it, i);
        _Centers.push_back(*it);
    }
}

int Kmeans::NearestCenter(const Point& p) {
    float minDistance = std::numeric_limits<float>::max();
    int k_id = -1;
    auto it = _Centers.begin();
    for (int k = 0; it != _Centers.end(); ++it, ++k) {
        float dist = Distance(p, *it);
        if (dist < minDistance) {
            minDistance = dist;
            k_id = k;
        }
    }
    return k_id;
}

void Kmeans::Cluster() {
#ifdef USING_OMP
#pragma omp parallel for
#endif
    for (auto& point : _Points) {
        point._group = NearestCenter(point);
    }
}

void Kmeans::Center() {
    std::vector<Point> centers(_K, Point(0.0, 0.0, 0.0)); 
    std::vector<int> counts(_K, 0);

    for (const auto& point : _Points) {
        for (int i = 0; i < 3; ++i) {
            centers[point._group][i] += point[i];
        }
        counts[point._group]++;
    }

    auto it = _Centers.begin();
    for (int i = 0; i < _K; ++i, ++it) {
        if (counts[i] > 0) {
            for (int j = 0; j < 3; ++j) {
                centers[i][j] /= counts[i];
            }
        }
        *it = centers[i];
    }
}

float Kmeans::CalculateTotalClusteringError() {
    float totalWithinClusterError = 0;
    float totalBetweenClusterDistance = 0;

    // Calculate within-cluster error
    for (const auto& point : _Points) {
        auto it = _Centers.begin();
        std::advance(it, point._group);  // Advance iterator to the correct center based on group index
        const Point& centroid = *it;

        // // Manually calculate the Euclidean distance
        // float dist = 0.0;
        // for (int i = 0; i < 3; ++i) {
        //     dist += pow(centroid[i] - point[i], 2);
        // }
        // dist = sqrt(dist); // Take the square root to get the Euclidean distance

        // totalWithinClusterError += dist;

        totalWithinClusterError += Kmeans::Distance(centroid, point);
        
    }

    // Calculate between-cluster distance
    for (auto it1 = _Centers.begin(); it1 != _Centers.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != _Centers.end(); ++it2) {
            totalBetweenClusterDistance += Distance(*it1, *it2);
        }
    }

    // Typically, you might subtract between-cluster distance to emphasize separation, 
    // but here we'll just add it to stick with the error definition provided
    return totalWithinClusterError - totalBetweenClusterDistance;

}

void Kmeans::RunKmean() {
    std::list<Point> oldCenters;
    _MaxIteration = 66;
    for (int iteration = 0; iteration < _MaxIteration; iteration++) {
        oldCenters = _Centers;
        Cluster();
        Center();
        float clusteringError = CalculateTotalClusteringError();
        std::cout << "Iteration " << iteration << " Clustering Error: " << clusteringError << std::endl;
        if (std::abs(clusteringError) < 0.0001) {
            break;
        }
    }
}

void Kmeans::PrintPointList(const std::list<Point>& pointList) {
    for (const auto& point : pointList) {
        std::cout << point[0] << "\t" << point[1] << "\t" << point[2] << std::endl;
    }
}









// int main() {
//     string filepath = "/Users/OctoPii/Desktop/GMDA_RPTREE/GLU_dihedrals.csv";
//     ifstream file(filepath);
//     list<Point> data;
//     string line;
    
//     int N_SAMPLES = 1000; 

//     cout << "File found successfully.";
//     cout << "STARTING READ AND PROCESSING DATA." << endl;

//     if (file.is_open()) {
//         // Skip the first line (header)
//         getline(file, line);

//         while (getline(file, line)) {
//             stringstream ss(line);
//             vector<double> values;
//             string cell;

//             // Read values for columns "phi", "psi", and "omega"
//             for (int i = 0; i < 3; ++i) {
//                 if (!getline(ss, cell, ',')) {
//                     cerr << "Error: Incomplete line." << endl;
//                     return 1;
//                 }

//                 try {
//                     values.push_back(stod(cell));
//                 } catch (const invalid_argument& e) {
//                     cerr << "Failed to convert '" << cell << "' to double: " << e.what() << endl;
//                     cerr << "Problematic line: " << line << endl;
//                     return 1;
//                 }
//             }

//             if (!values.empty()) {
//                 Point p(values[0], values[1], values[2]);  // Create a Point for each line
//                 data.push_back(p);
//             }
//             if (data.size() == N_SAMPLES) break;
//         } 
//         file.close();
//     } else {
//         cout << "Failed to open file: " << filepath << endl;
//         return 1;
//     }

//     // Display the first five points from data
//     auto it = data.begin();
//     int count = 0;
//     for (; it != data.end() && count < 5; ++it, ++count) {
//         cout << "Phi: " << (*it)[0] << ", Psi: " << (*it)[1] << ", Omega: " << (*it)[2] << endl;
//     }

//     int MINSIZE = 10; 
//     int DEPTH_MAX = int((log(N_SAMPLES/MINSIZE))/log(2))+1;
//     int N_CLUSTER_MAX = pow(DEPTH_MAX, 2);
//     int k = N_CLUSTER_MAX; // Number of clusters
//     Kmeans kmeans(k, data.size());
//     kmeans.InitPoints(data); // Initialize Kmeans with the data
//     kmeans.InitCenters(); // Initialize centers

//     std::cout << "Kmeans initialized. Computing error ..." << std::endl;

//     // Timer
//     auto start = std::chrono::high_resolution_clock::now();
//     kmeans.RunKmean();
//     auto stop = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<chrono::milliseconds>(stop - start);
    
//     cout << "Total clustering error: " << kmeans.CalculateTotalClusteringError() << endl;

//     double durationInMinutes = duration.count() / 60000.0;  // Convert to minutes
//     cout << "Time taken: " << std::fixed << std::setprecision(3) << durationInMinutes << " min." << endl;

//     return 0;
// }