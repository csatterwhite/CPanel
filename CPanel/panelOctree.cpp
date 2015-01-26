//
//  panelOctree.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panelOctree.h"

Eigen::Vector3d panelOctree::findRefPoint(const panel &obj)
{
    return obj.getCenter();
}