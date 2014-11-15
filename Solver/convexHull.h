//
//  convexHull.h
//  CPanel
//
//  Created by Chris Satterwhite on 9/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__convexHull__
#define __CPanel__convexHull__

#include <stdio.h>
#include <Eigen/Dense>
#include <vector>
#include <list>
#include <algorithm>
#include "math.h"


class convexHull
{
    struct member
    {
        member(const Eigen::Vector3d &point, const Eigen::Vector3d &basePnt);
        double x,y,z,theta,d;
    };
    struct compareTheta
    {
        inline bool operator()(member* p1, member* p2)
        {
            return (p1->theta < p2->theta);
        }
    };
    struct compareDistAscending
    {
        inline bool operator()(member* p1, member* p2)
        {
            return (p1->d < p2->d);
        }
    };
    struct compareDistDescending
    {
        inline bool operator()(member* p1, member* p2)
        {
            return (p1->d > p2->d);
        }
    };
    
    std::vector<member*> members;
    std::vector<member*> hull;
    
    void computeHull();
    bool boundary;
    // TRUE : Points on boundary consider inside hull and therefore not included in hull vector.
    // FALSE : Points on boundary considered outside of hull and therefore included in hull vector
    
    Eigen::Vector3d makeVector(member* p1, member* p2);
    
public:
    
    convexHull(Eigen::MatrixXd points, bool boundary);
    ~convexHull()
    {
        for (int i = 0; i<members.size(); i++)
        {
            delete members[i];
        }
    }
    std::vector<member*> getHull() {return hull;}
};

#endif /* defined(__CPanel__convexHull__) */
