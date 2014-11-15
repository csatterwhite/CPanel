//
//  convexHullTests.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 9/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "convexHullTests.h"

convexHullTests::convexHullTests()
{
    TEST_ADD(convexHullTests::testPntInside);
    TEST_ADD(convexHullTests::testPntOutside);
    TEST_ADD(convexHullTests::testPntOnBoundary);
}

void convexHullTests::testPntInside()
{
    Eigen::Vector3d p1;
    p1 << 0,1,0;
    Eigen::Vector3d p2;
    p2 << 1,1,0;
    Eigen::Vector3d p3;
    p3 << 0,0,0;
    Eigen::Vector3d p4;
    p4 << 0.333,0.666,0;
    
    Eigen::MatrixXd points(4,3);
    points.row(0) = p1;
    points.row(1) = p2;
    points.row(2) = p3;
    points.row(3) = p4;
    
    convexHull testHull(points,true);
    
    TEST_ASSERT(testHull.getHull().size() == 3);
}

void convexHullTests::testPntOutside()
{
    Eigen::Vector3d p1;
    p1 << 0,1,0;
    Eigen::Vector3d p2;
    p2 << 1,1,0;
    Eigen::Vector3d p3;
    p3 << 0,0,0;
    Eigen::Vector3d p4;
    p4 << 1,0,0;
    
    Eigen::MatrixXd points(4,3);
    points.row(0) = p1;
    points.row(1) = p2;
    points.row(2) = p3;
    points.row(3) = p4;
    
    convexHull testHull(points,true);
    
    TEST_ASSERT(testHull.getHull().size() == 4);
}

void convexHullTests::testPntOnBoundary()
{
    Eigen::Vector3d p1;
    p1 << 0,0.5,0;
    Eigen::Vector3d p2;
    p2 << 1,1,0;
    Eigen::Vector3d p3;
    p3 << 0,0,0;
    Eigen::Vector3d p4;
    p4 << 0.5,0.5,0;
    Eigen::Vector3d p5;
    p5 << 0,1,0;
    Eigen::Vector3d p6;
    p6 << 0.5,1,0;
    
    Eigen::MatrixXd points(6,3);
    points.row(0) = p1;
    points.row(1) = p2;
    points.row(2) = p3;
    points.row(3) = p4;
    points.row(4) = p5;
    points.row(5) = p6;
    
    convexHull testBoundaryInside(points,true);
    convexHull testBoundaryOutside(points,false);
    
    TEST_ASSERT(testBoundaryInside.getHull().size() == 3);
    TEST_ASSERT(testBoundaryOutside.getHull().size() == 6);
}

