//
//  PanelTest.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 3/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "PanelTest.h"


PanelTest::PanelTest()
{
    TEST_ADD(PanelTest::test_setPnts);
    TEST_ADD(PanelTest::test_addNeighbor);
}


void PanelTest::test_setPnts()
{
    panel p;
    Eigen::Vector3i vertices = {1,2,3};
    Eigen::MatrixXd nodes(3,3);
    nodes << 0,0,0,3,0,0,0,3,0;
    p.setGeom(vertices,nodes);
    
    bool flag = true;
    Eigen::Vector3i verts = p.getVerts();
    if (verts(0) != 0 || verts(1) != 1 || verts(2) != 2)
    {
        flag = false;
    }
    TEST_ASSERT_MSG(flag, "Vertices Set Incorrectly")
    
    flag = true;
    Eigen::Vector3d center = p.getCenter();
    if (center(0) != 1 || center(1) != 1 || center(2) != 0)
    {
        flag = false;
    }
    TEST_ASSERT_MSG(flag, "Center Set Incorrectly");
    
    flag = true;
    Eigen::Vector3d norm = p.getNormal();
    if (norm(0) != 0 || norm(1) != 0 || norm(2) != 1)
    {
        flag = false;
    }
    TEST_ASSERT_MSG(flag,"Normal Set Incorrectly");
}

