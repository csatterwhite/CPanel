//
//  horseshoeVortex.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/10/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "horseshoeVortex.h"

double horseshoeVortex::phiInfluence(Eigen::Vector3d POI)
{
    double x,xa,y,ya,yb,z,term1,term2,term3,term4;
    x = POI(0);
    xa = (p1(0)+p2(0))/2;
    y = POI(1);
    ya = p1(1);
    yb = p2(1);
    
    Eigen::Vector3d xdir;
    xdir << 1,0,0;
    Eigen::Vector3d t = (p2-p1).normalized();
    Eigen::Vector3d normal = xdir.cross(t);
    z = (POI-((p2-p1)/2)).dot(normal);
    
    
    term1 = atan(z/(y-yb));
    term2 = -atan(z/(y-ya));
    term3 = atan((yb-y)*(x-xa)/(z*sqrt(pow(x-xa,2)+ pow(y-yb,2)+pow(z,2))));
    term4 = -atan((ya-y)*(x-xa)/(z*sqrt(pow(x-xa,2)+pow(y-ya,2)+pow(z,2))));
    double influence = -1/(4*M_PI)*(term1+term2+term3+term4);
    return influence;
}

double horseshoeVortex::horseshoePhi(Eigen::Vector3d POI)
{
    return gamma*phiInfluence(POI);
}