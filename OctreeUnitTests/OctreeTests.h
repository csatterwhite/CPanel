//
//  OctreeTests.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__OctreeTests__
#define __CPanel__OctreeTests__

#include <iostream>
#include "cpptest/cpptest.h"
#include "octree.h"
#include "node.h"
#include "testOctreeClass.h"

class OctreeTests : public Test::Suite
{
    std::vector<testObj*> data;
    
    void test_constructor();
    void test_setMaxMembers();
    void test_addData();
//    void test_findNodeContainingMember();
    
public:
    OctreeTests();
    
    ~OctreeTests()
    {
        for (int i=0; i<data.size(); i++)
        {
            delete data[i];
        }
    }
    
};

#endif /* defined(__CPanel__OctreeTests__) */
