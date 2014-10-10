//
//  InfluenceTests.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/7/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "InfluenceTests.h"


influenceTests::influenceTests()
{
    TEST_ADD(testPntClose)
    TEST_ADD(testPntFar)
}

influenceTests::testPntClose()
{
    panel p(1);
    Eigen::Matrix3d nodes;
    double a = 1;
    double h = a*sin(M_PI/3);
    nodes(0,0) = 2*h/3;
    nodes(0,1) = 0;
    nodes(0,2) = 0;
    nodes(1,0) = -h/3;
    nodes(1,1) = a/2;
    nodes(1,2) = 0;
    nodes(2,0) = -h/3;
    nodes(2,1) = -a/2;
    nodes(2,2) = 0;
    p.setGeom(verts,nodes);
    
    double phiSource;
    Eigen::Vector3d vSource;
    Eigen::Vector3d POI;
    POI << 2,1,1;
    p.sourceInfluence(1,POI,nodes,phiSource,vSource);
    
    
}