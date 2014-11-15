//
//  PanelTest.h
//  CPanel
//
//  Created by Chris Satterwhite on 3/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__PanelTest__
#define __CPanel__PanelTest__

#include <iostream>
#include <cpptest/cpptest.h>
#include "panel.h"
#include "panelOctree.h"


class PanelTest : public Test::Suite
{
    void test_setGeomTri();
    void test_setGeomQuad();
    void test_addNeighbor();

public:
    PanelTest();
    
};

#endif /* defined(__CPanel__PanelTest__) */
