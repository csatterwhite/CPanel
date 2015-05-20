//
//  convexHull.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 9/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "convexHull.h"

convexHull::convexHull(Eigen::MatrixXd points, bool bound) : boundary(bound)
{
    Eigen::Vector3d basePnt = points.row(0);
    double minTheta = M_PI;
    double maxTheta = 0;
    
    for (int i=1; i<points.rows(); i++)
    {
        if ((points(i,1) == basePnt(1) && points(i,0) < basePnt(0)) || (points(i,1) < basePnt(1)))
        {
            basePnt = points.row(i);
        }
    }
    for (int i=0; i<points.rows(); i++)
    {
        member* m = new member(points.row(i),basePnt);
        if (points(i,0) != basePnt(0) || points(i,1) != basePnt(1) || points(i,2) != basePnt(2))
        {
            if (m->theta < minTheta)
            {
                minTheta = m->theta;
            }
            if (m->theta > maxTheta)
            {
                maxTheta = m->theta;
            }
        }
        if (points(i,0) == basePnt(0) && points(i,1) == basePnt(1) && points(i,2) == basePnt(2))
        {
            members.insert(members.begin(),m);
        }
        else
        {
            members.push_back(m);
        }
    }
    std::sort(members.begin()+1,members.end(),compareTheta());
    int minBegin = -1;
    int minEnd = -1;
    unsigned long maxBegin = members.size();
    for (int i=1; i<members.size()-1; i++)
    {
        if (members[i]->theta == minTheta && minBegin == -1)
        {
            minBegin = i;
        }
        else if (members[i]->theta != minTheta && minBegin != -1 && minEnd == -1)
        {
            minEnd = i;
        }
        else if (members[i]->theta == maxTheta && maxBegin == members.size())
        {
            maxBegin = i;
        }
    }
    std::sort(members.begin()+minBegin,members.begin()+minEnd,compareDistAscending());
    std::sort(members.begin()+maxBegin,members.end(),compareDistDescending());
    computeHull();
}


convexHull::member::member(const Eigen::Vector3d &point, const Eigen::Vector3d &basePnt)
{
    x = point(0);
    y = point(1);
    z = point(2);
    theta = atan2((y-basePnt(1)),(x-basePnt(0)));
    d = sqrt(pow(y-basePnt(1),2)+pow(x-basePnt(0),2));
}

void convexHull::computeHull()
{
    hull.push_back(members.front());
    for (int i=0; i<members.size()-1; i++)
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

bool convexHull::compareNodes(std::vector<Eigen::Vector3d> nodesLocal)
{
    if (hull.size() > nodesLocal.size())
    {
        return false;
    }
    
    Eigen::Vector3d pMember;
    bool breakFlag = false;
    
    for (int i=0; i<hull.size(); i++)
    {
        pMember(0) = hull[i]->x;
        pMember(1) = hull[i]->y;
        pMember(2) = hull[i]->z;
        for (int j=0; j<nodesLocal.size(); j++)
        {
            if (pMember == nodesLocal[j])
            {
                breakFlag = true;
                break;
            }
        }
        if (breakFlag)
        {
            breakFlag = false;
            continue;
        }
        else
        {
            return false;
        }
    }
    return true;
}