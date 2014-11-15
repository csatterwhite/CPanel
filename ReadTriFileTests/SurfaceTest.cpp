//
//  SurfaceTest.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 3/27/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "SurfaceTest.h"


SurfaceTest::SurfaceTest()
{
    TEST_ADD(SurfaceTest::test_addPanel);
}

void SurfaceTest::test_addPanel()
{
    Eigen::Vector3i indices = {0,1,2};
    Eigen::MatrixXd nodes(3,3);
    nodes << 0,0,0,0,1,0,1,0,0;
    surface surf(1,&nodes);
    surf.addPanel(indices);
    TEST_ASSERT_MSG(surf.getPanels().size() == 1,"Failed to add panel to surface");
    TEST_ASSERT_MSG(surf.getPanels()[0] != NULL, "Panel not created (pointer is NULL)");
}

