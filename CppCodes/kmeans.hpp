#ifndef _KMEANS_H
#define _KMEANS_H

#include<iostream> 
#include<list>
#include<vector>
#include<cmath>
#include<fstream>
#include<limits>
#include<map>
#include<random>

#ifdef USING_OMP
#include <omp.h>
#endif

struct Point {
    float coords[3]; // Store x, y, z in an array
    int _group;

    Point(float x = 0, float y = 0, float z = 0, int g = 0) : _group(g) {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }

    float& operator[](int index) { return coords[index]; }
    const float& operator[](int index) const { return coords[index]; }
};

class Kmeans {
public:
    Kmeans(int k, int pointnumber = 0);
    void InitPoints();
    void InitPoints(const std::list<Point>& pointlist);
    void InitPoints(std::string fileName);
    virtual void InitCenters();
    void InitSpecifiedCenters();
    int NearestCenter(const Point& p);
    float CalculateTotalClusteringError();
    void RunKmean();
    void Cluster();
    void Center();
    float Distance(const Point& p1, const Point& p2);
    void PrintPointList(const std::list<Point>& pointList);

// private:
    std::list<Point> _Points;
    std::list<Point> _Centers;
    int _MaxIteration;
    int _K;
    int _PointNumber;

};

#endif