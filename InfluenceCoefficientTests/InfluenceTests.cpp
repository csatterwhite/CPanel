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
    TEST_ADD(influenceTests::testPntClose)
    TEST_ADD(influenceTests::testPntFar)
    TEST_ADD(influenceTests::testCollocationPnt)
}

void influenceTests::testPntClose()
{
    Eigen::Vector3d POI;
    POI << 2,1,1; // |POI|/(side length) ~ 2.5 so panel formulation will be used;
    
    panel p(1);
    double a = 1;
    double h = a*sin(M_PI/3);
    Eigen::Matrix3d nodes;
    nodes(0,0) = 2*h/3;
    nodes(0,1) = 0;
    nodes(0,2) = 0;
    nodes(1,0) = -h/3;
    nodes(1,1) = a/2;
    nodes(1,2) = 0;
    nodes(2,0) = -h/3;
    nodes(2,1) = -a/2;
    nodes(2,2) = 0;
    Eigen::Vector3i verts;
    verts << 0,1,2;
    p.setGeom(verts,nodes);
    
    double phiSource;
    Eigen::Vector3d vSource;
    p.sourceInfluence(1,POI,nodes,phiSource,vSource);
    TEST_ASSERT(phiSource > -0.01409318 && phiSource < -0.01409317)
    TEST_ASSERT(vSource(0) > 0.00469304 && vSource(0) < 0.00469305)
    TEST_ASSERT(vSource(1) > 0.00235937 && vSource(1) < 0.00235938)
    TEST_ASSERT(vSource(2) > 0.00240016 && vSource(2) < 0.00240017)
    
    double phiDoublet;
    Eigen::Vector3d vDoublet;
    p.doubletInfluence(1, POI, nodes, phiDoublet, vDoublet);
    TEST_ASSERT(phiDoublet > -0.00240017 && phiDoublet < -0.00240016)
    TEST_ASSERT(vDoublet(0) > 0.00241728 && vDoublet(0) < 0.00241729)
    TEST_ASSERT(vDoublet(1) > 0.00122408 && vDoublet(1) < 0.00122409)
    TEST_ASSERT(vDoublet(2) > -0.00114367 && vDoublet(2) < -0.00114366)
}

void influenceTests::testPntFar()
{
    Eigen::Vector3d POI;
    POI << 4,3,2; // |POI|/(side length) ~ 5.4 so far field approximation will be used;
    
    panel p(1);
    double a = 1;
    double h = a*sin(M_PI/3);
    Eigen::Matrix3d nodes;
    nodes(0,0) = 2*h/3;
    nodes(0,1) = 0;
    nodes(0,2) = 0;
    nodes(1,0) = -h/3;
    nodes(1,1) = a/2;
    nodes(1,2) = 0;
    nodes(2,0) = -h/3;
    nodes(2,1) = -a/2;
    nodes(2,2) = 0;
    Eigen::Vector3i verts;
    verts << 0,1,2;
    p.setGeom(verts,nodes);
    
    double phiSource;
    Eigen::Vector3d vSource;
    p.sourceInfluence(1,POI,nodes,phiSource,vSource);
    TEST_ASSERT(phiSource > -0.00639871 && phiSource < -0.00639870)
    TEST_ASSERT(vSource(0) > 0.00088257 && vSource(0) < 0.00088258)
    TEST_ASSERT(vSource(1) > 0.00066193 && vSource(1) < 0.00066194)
    TEST_ASSERT(vSource(2) > 0.00044128 && vSource(2) < 0.00044129)
    
    double phiDoublet;
    Eigen::Vector3d vDoublet;
    p.doubletInfluence(1, POI, nodes, phiDoublet, vDoublet);
    TEST_ASSERT(phiDoublet > -0.00044129 && phiDoublet < -0.00044128)
    TEST_ASSERT(vDoublet(0) > 0.00018260 && vDoublet(0) < 0.00018261)
    TEST_ASSERT(vDoublet(1) > 0.00013695 && vDoublet(1) < 0.00013696)
    TEST_ASSERT(vDoublet(2) > -0.00012935 && vDoublet(2) < -0.00012934)
}

void influenceTests::testCollocationPnt()
{
    Eigen::Vector3d POI;
    POI << 0,0,0;
    
    panel p(1);
    double a = 1;
    double h = a*sin(M_PI/3);
    Eigen::Matrix3d nodes;
    nodes(0,0) = 2*h/3;
    nodes(0,1) = 0;
    nodes(0,2) = 0;
    nodes(1,0) = -h/3;
    nodes(1,1) = a/2;
    nodes(1,2) = 0;
    nodes(2,0) = -h/3;
    nodes(2,1) = -a/2;
    nodes(2,2) = 0;
    Eigen::Vector3i verts;
    verts << 0,1,2;
    p.setGeom(verts,nodes);
    
    double phiSource;
    Eigen::Vector3d vSource;
    p.sourceInfluence(1,POI,nodes,phiSource,vSource);
    TEST_ASSERT(phiSource == 0)
    TEST_ASSERT(vSource(0) == 0)
    TEST_ASSERT(vSource(1) == 0)
    TEST_ASSERT(vSource(2) == 0.5)
    
    double phiDoublet;
    Eigen::Vector3d vDoublet;
    p.doubletInfluence(1,POI,nodes,phiDoublet,vDoublet);
    TEST_ASSERT(phiDoublet == -0.5)
    TEST_ASSERT(vDoublet(0) == 0)
    TEST_ASSERT(vDoublet(1) == 0)
    TEST_ASSERT(vDoublet(2) > 1.432394487 && vDoublet(2) < 1.432394488)
    
    
}