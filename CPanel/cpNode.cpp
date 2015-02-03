//
//  cpNode.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 2/2/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "cpNode.h"
#include "edge.h"


cpNode::cpNode(Eigen::Vector3d pnt,int index) : pnt(pnt), index(index) {}

Eigen::Vector3d cpNode::operator-=(const cpNode &rhs)
{
    return pnt-rhs.getPnt();
}

Eigen::Vector3d operator-(cpNode lhs, const cpNode &rhs)
{
    return lhs -= rhs;
}

Eigen::Vector3d cpNode::operator+=(const cpNode &rhs)
{
    return pnt+rhs.getPnt();
}

Eigen::Vector3d operator+(cpNode lhs, const cpNode &rhs)
{
    return lhs += rhs;
}

void cpNode::addEdge(edge* e)  {edges.push_back(e);}