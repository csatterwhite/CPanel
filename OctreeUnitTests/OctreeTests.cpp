//
//  OctreeTests.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 4/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "OctreeTests.h"

OctreeTests::OctreeTests()
{
    TEST_ADD(OctreeTests::test_constructor);
    TEST_ADD(OctreeTests::test_setMaxMembers);
    TEST_ADD(OctreeTests::test_addData);
//    TEST_ADD(OctreeTests::test_findNodeContainingMember);
}

void OctreeTests::test_constructor()
{
    testOctreeClass testOctree;
    TEST_ASSERT(testOctree.getMaxMembersPerNode()==10)
}

void OctreeTests::test_setMaxMembers()
{
    int max = 5;
    testOctreeClass testOctree;
    testOctree.setMaxMembers(max);
    TEST_ASSERT(testOctree.getMaxMembersPerNode()==max)
}

void OctreeTests::test_addData()
{
    int nX = 3;
    int nY = 3;
    int nZ = 3;
    testObj* obj;
    for (int i=0; i<nX; i++)
    {
        for (int j=0; j<nY; j++)
        {
            for (int k=0; k<nZ; k++)
            {
                std::array<double,3> center;
                center[0] = i;
                center[1] = j;
                center[2] = k;
                obj = new testObj(center);
                data.push_back(obj);
            }
        }
    }
    testOctreeClass testOctree;
    testOctree.addData(data);
    
    TEST_ASSERT(testOctree.getMembers().size()==nX*nY*nZ);
    
    bool flag = true;
    for (int i=0; i<3; i++)
    {
        if (testOctree.getRootNode()->getOrigin()[i] != 1)
        {
            flag = false;
        }
    }
    TEST_ASSERT(flag)
    TEST_ASSERT(testOctree.getRootNode()->getParent() == NULL);
    
    node<testObj>* oldRoot = testOctree.getRootNode();
    
    testObj* newData;
    std::array<double,3> newPoint = {-2,-2,-2};
    newData = new testObj(newPoint);
    data.push_back(newData);
    testOctree.addData(newData);
    
    TEST_ASSERT(testOctree.getMembers().size()==nX*nY*nZ+1);
    TEST_ASSERT(testOctree.getRootNode() != oldRoot);
    
}