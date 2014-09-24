//
//  convexHull.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 9/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "convexHull.h"

convexHull::convexHull(std::vector<Eigen::Vector3d> points, bool boundary) : boundary(boundary)
{
    Eigen::Vector3d basePnt = points[0];
    for (int i=1; i<points.size(); i++)
    {
        if ((points[i](1) == basePnt(1) && points[i](0) < basePnt(0)) || (points[i](1) < basePnt(1)))
        {
            basePnt = points[i];
        }
    }
    for (int i=0; i<points.size(); i++)
    {
        member* m = new member(points[i],basePnt);
        if (points[i] == basePnt)
        {
            members.insert(members.begin(),m);
        }
        else
        {
            members.push_back(m);
        }
    }
    std::sort(members.begin()+1,members.end(),compareD());
    std::sort(members.begin()+1,members.end(),compareTheta());
}

convexHull::member::member(Eigen::Vector3d point, Eigen::Vector3d basePnt)
{
    x = point(0);
    y = point(1);
    y = point(2);
    theta = atan2((y-basePnt(1)),(x-basePnt(0)));
    d = sqrt(pow(y-basePnt(1),2)+pow(x-basePnt(0),2));
}

void convexHull::computeHull()
{
    hull.push_back(members.front());
    for (int i=0; i<members.size(); i++)
    {
        Eigen::Vector3d v1,v2;
        if (i==members.size()-2)
        {
            v1 = makeVector(members[i],members[i+1]);
            v2 = makeVector(members[i],members[0]);
        }
        else
        {
            v1 = makeVector(members[i],members[i+1]);
            v2 = makeVector(members[i],members[i+2]);
        }
        
        double crossZ = v1.cross(v2)(2);
        
        if ((crossZ > 0) || (crossZ == 0 && !boundary))
        {
            hull.push_back(members[i+1]);
        }
    }
}

Eigen::Vector3d convexHull::makeVector(member* p1, member* p2)
{
    Eigen::Vector3d vec;
    vec(0) = p2->x-p1->x;
    vec(1) = p2->y-p1->y;
    vec(2) = 0; //Only computing 2D complex hull
    return vec;
}