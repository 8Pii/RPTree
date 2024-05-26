#ifndef KMEANS_PLUS_PLUS_H
#define KMEANS_PLUS_PLUS_H

#include "kmeans.hpp"  // Make sure the include path is correct

class KmeansPlusPlus : public Kmeans {
public:
    KmeansPlusPlus(int k, int pointnumber);
    KmeansPlusPlus(int k);

    void InitCenters() override;  // Ensure override is correctly specified
    int NearestCenter(Point &p, int alreadyInitCenterNumber, float &minDistance);
};

#endif
