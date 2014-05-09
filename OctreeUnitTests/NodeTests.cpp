//
//  NodeTests.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 4/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "NodeTests.h"


NodeTests::NodeTests() : testNode(nullptr,std::array<double,3> {1,1,1},std::array<double,3> {1,1,1},0,10)
{
    TEST_ADD(NodeTests::test_constructor)
    TEST_ADD(NodeTests::test_isLeafNode);
    TEST_ADD(NodeTests::test_addMember);
    TEST_ADD(NodeTests::test_getNodeContainingMember);
    TEST_ADD(NodeTests::test_createParent);
}

void NodeTests::test_constructor()
{
    // Test origin is set correctly
    bool flag = true;
    for (int i=0; i<3; i++)
    {
        if (testNode.getOrigin()[i] != 1)
        {
            flag = false;
            break;
        }
    }
    TEST_ASSERT(flag);
    
    // Test halfDimension is set correctly
    flag = true;
    for (int i=0; i<3; i++)
    {
        if (testNode.getHalfDimension()[i] != 1)
        {
            flag = false;
            break;
        }
    }
    TEST_ASSERT(flag);
    
    // Test children are NULL
    flag = true;
    for (int i=0; i<8; i++)
    {
        if (testNode.getChild(i) != nullptr)
        {
            flag = false;
        }
    }
    TEST_ASSERT(flag);
    
    // Test level set
    TEST_ASSERT(testNode.getLevel() == 1);
}

void NodeTests::test_isLeafNode()
{
    TEST_ASSERT(testNode.isLeafNode());
}

void NodeTests::test_addMember()
{
    int maxMembers = 10;
    testNode.setMaxMembers(maxMembers);
    
    int nX = 3;
    int nY = 3;
    int nZ = 3;
    int counter = 0;
    for (int i=0; i<nX; i++)
    {
        for (int j=0; j<nY; j++)
        {
            for (int k=0; k<nZ; k++)
            {
                point obj;
                obj[0] = i;
                obj[1] = j;
                obj[2] = k;
                member<point> member(&obj,obj);
                testNode.addMember(member);
                counter += 1;
                if (counter <= maxMembers)
                {
                    // Number of members in node should be the same as the number added.
                    TEST_ASSERT(testNode.getMembers().size() == counter);
                }
                else
                {
                    // Once maxMembers is exceeded, all members should be pushed to children and therefore node will no longer be leaf and should contain no members.
                    TEST_ASSERT(!testNode.isLeafNode());
                    TEST_ASSERT(testNode.getMembers().size() == 0);
                }
            }
        }
    }
}

void NodeTests::test_getNodeContainingMember()
{
    std::array<double,3> origin = testNode.getOrigin();
    std::array<double,3> obj;
    obj[0] = origin[0]-1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]-1;
    member<point> member0(&obj,obj);
    testNode.addMember(member0);
    TEST_ASSERT(testNode.getChildContainingMember(member0) == 0)
    
    obj[0] = origin[0]-1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]+1;
    member<point> member1(&obj,obj);
    testNode.addMember(member1);
    TEST_ASSERT(testNode.getChildContainingMember(member1) == 1)
    
    obj[0] = origin[0]-1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]-1;
    member<point> member2(&obj,obj);
    testNode.addMember(member2);
    TEST_ASSERT(testNode.getChildContainingMember(member2) == 2)
    
    obj[0] = origin[0]-1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]+1;
    member<point> member3(&obj,obj);
    testNode.addMember(member3);
    TEST_ASSERT(testNode.getChildContainingMember(member3) == 3)
    
    obj[0] = origin[0]+1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]-1;
    member<point> member4(&obj,obj);
    testNode.addMember(member4);
    TEST_ASSERT(testNode.getChildContainingMember(member4) == 4)
    
    obj[0] = origin[0]+1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]+1;
    member<point> member5(&obj,obj);
    testNode.addMember(member5);
    TEST_ASSERT(testNode.getChildContainingMember(member5) == 5)
    
    obj[0] = origin[0]+1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]-1;
    member<point> member6(&obj,obj);
    testNode.addMember(member6);
    TEST_ASSERT(testNode.getChildContainingMember(member6) == 6)
    
    obj[0] = origin[0]+1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]+1;
    member<point> member7(&obj,obj);
    testNode.addMember(member7);
    TEST_ASSERT(testNode.getChildContainingMember(member7) == 7)
}

void NodeTests::test_createParent()
{
    testNode.createParent(0);
    TEST_ASSERT(testNode.getChild(7) != nullptr);
    TEST_ASSERT(testNode.getLevel() == 2);
}

