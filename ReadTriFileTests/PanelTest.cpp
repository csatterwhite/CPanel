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
    panel p(1);
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

void PanelTest::test_addNeighbor()
{
    Eigen::MatrixXd nodes(6,3);
    Eigen::MatrixXi indices(4,3);
    std::vector<panel*> panels;
    
    nodes << 0,0,0 , 1,0,0 , 2,0,0 , 2,1,0 , 1,1,0 , 0,1,0;
    
    indices << 1,2,6 , 6,2,5 , 5,2,3 , 3,4,5;
    
    //  6 ___5___4
    //   |\  |\  |
    //   | \ | \ |
    //   |__\|__\|
    //  1    2    3
    
    panel* p;
    for (int i=0; i<4; i++)
    {
        p = new panel(1);
        p->setGeom(indices.row(i),nodes);
        panels.push_back(p);
    }
    
    panels[0]->checkNeighbor(panels[1]);
    panels[0]->checkNeighbor(panels[2]);
    panels[0]->checkNeighbor(panels[3]);
    
    TEST_ASSERT_MSG(panels[0]->getNeighbors().size() == 1, "Wrong number of primary neighbors")
    TEST_ASSERT_MSG(panels[1]->getNeighbors().size() == 1, "Reciprocal addition of neighbor failed")
    TEST_ASSERT_MSG(panels[2]->getNeighbors().size() == 0, "False neighbor added to panel")
    TEST_ASSERT_MSG(panels[3]->getNeighbors().size() == 0, "False neighbor added to panel")
    
    panels[1]->checkNeighbor(panels[0]);
    
    TEST_ASSERT_MSG(panels[1]->getNeighbors().size() == 1, "Duplicate neighbor added")
    
    
    for (int i=0; i<4; i++)
    {
        delete panels[i];
    }
}

