//
//  testOctreeClass.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__testOctreeClass__
#define __CPanel__testOctreeClass__

#include <iostream>
#include "octree.h"
#include <array>
#include "testObj.h"

class testOctreeClass : public octree<testObj>
{
public:
    std::array<double,3> findRefPoint(const testObj* member) {return member->getCenter();}
    
    testOctreeClass() : octree() {}
    
};

#endif /* defined(__CPanel__testOctreeClass__) */
