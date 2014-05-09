//
//  NodeTests.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__NodeTests__
#define __CPanel__NodeTests__

#include <iostream>
#include "cpptest/cpptest.h"
#include "node.h"
#include "member.h"

class NodeTests : public Test::Suite
{
    typedef std::array<double,3> point;
    
    node<std::array<double,3>> testNode;
    void test_constructor();
    void test_isLeafNode();
    void test_addMember();
    void test_getNodeContainingMember();
    void test_createParent();

public:
    NodeTests();
};

#endif /* defined(__CPanel__NodeTests__) */
