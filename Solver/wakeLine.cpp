//
//  wakeLine.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakeLine.h"

wakeLine::wakeLine(bodyPanel* upper, bodyPanel* lower,Eigen::Vector3d normal) : upper(upper), lower(lower), normal(normal)
{
    setDimensions();
}

double wakeLine::getStrength()
{
    return upper->getMu()-lower->getMu();
}

void wakeLine::setDimensions()
{
    std::vector<Eigen::Vector3d> TEverts;
    Eigen::MatrixXd* nodes = upper->getNodes();
    Eigen::VectorXi uVerts = upper->getVerts();
    Eigen::VectorXi lVerts = lower->getVerts();
    for (int i=0; i<upper->getVerts().size(); i++)
    {
        for (int j=0; j<lower->getVerts().size(); j++)
        {
            if (nodes->row(uVerts(i)) == nodes->row(lVerts(j)))
            {
                TEverts.push_back(nodes->row(uVerts(i)));
            }
        }
    }
    std::sort(TEverts.begin(), TEverts.end(), [](Eigen::Vector3d p1, Eigen::Vector3d p2) {return p1(1)<p2(1);});
    p1 = TEverts[0];
    p2 = TEverts[1];
    pMid = (p1+p2)/2;
}
