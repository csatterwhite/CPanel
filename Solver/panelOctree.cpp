//
//  panelOctree.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panelOctree.h"

std::array<double,3> panelOctree::findRefPoint(const panel &obj)
{
    std::array<double,3> center;
    for (int i=0; i<3; i++)
    {
        center[i] = obj.getCenter()(i);
    }
    return center;
}