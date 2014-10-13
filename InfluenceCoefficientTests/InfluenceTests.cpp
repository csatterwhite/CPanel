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
    TEST_ASSERT(isEqual(phiSource,-0.014093178,9))
    TEST_ASSERT(isEqual(vSource(0),0.004693044,9))
    TEST_ASSERT(isEqual(vSource(1),0.002359374,9))
    TEST_ASSERT(isEqual(vSource(2),0.002400169,9))
    
    double phiDoublet;
    Eigen::Vector3d vDoublet;
    p.doubletInfluence(1, POI, nodes, phiDoublet, vDoublet);
    TEST_ASSERT(isEqual(phiDoublet,-0.002400169,9))
    TEST_ASSERT(isEqual(vDoublet(0),0.002417289,9))
    TEST_ASSERT(isEqual(vDoublet(1),0.001224085,9))
    TEST_ASSERT(isEqual(vDoublet(2),-0.001143661,9))
    
    double IC;
    bool dirichletBC = true;
    IC = p.sourceIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,phiSource,9))
    IC = p.doubletIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,phiDoublet,9));
    
    dirichletBC = false;
    IC = p.sourceIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,vSource(2),9))
    IC = p.doubletIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,vDoublet(2),9));
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
    TEST_ASSERT(isEqual(phiSource,-0.006398700,9))
    TEST_ASSERT(isEqual(vSource(0),0.000882579,9))
    TEST_ASSERT(isEqual(vSource(1),0.000661934,9))
    TEST_ASSERT(isEqual(vSource(2),0.000441289,9))
    
    double phiDoublet;
    Eigen::Vector3d vDoublet;
    p.doubletInfluence(1, POI, nodes, phiDoublet, vDoublet);
    TEST_ASSERT(isEqual(phiDoublet,-0.000441289,9))
    TEST_ASSERT(isEqual(vDoublet(0),0.000182602,9))
    TEST_ASSERT(isEqual(vDoublet(1),0.000136952,9))
    TEST_ASSERT(isEqual(vDoublet(2),-0.000129343,9))
    
    double IC;
    bool dirichletBC = true;
    IC = p.sourceIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,phiSource,9))
    IC = p.doubletIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,phiDoublet,9));
    
    dirichletBC = false;
    IC = p.sourceIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,vSource(2),9))
    IC = p.doubletIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,vDoublet(2),9));
    
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
    TEST_ASSERT(isEqual(vDoublet(2),1.432394487,9))
    
    double IC;
    bool dirichletBC = true;
    IC = p.sourceIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,phiSource,9))
    IC = p.doubletIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,phiDoublet,9));
    
    dirichletBC = false;
    IC = p.sourceIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,vSource(2),9))
    IC = p.doubletIC(POI,nodes,dirichletBC);
    TEST_ASSERT(isEqual(IC,vDoublet(2),9));
    
    
}

bool influenceTests::isEqual(double var1, double var2, int decimalPrecision)
{
    double eps = pow(10,-decimalPrecision);
    return (abs(var1-var2)<eps);
    
}
