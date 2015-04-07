//
//  cpNode.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 2/2/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "cpNode.h"
#include "edge.h"


cpNode::cpNode(Eigen::Vector3d pnt,int index) : pnt(pnt), index(index), TEnode(false) {}

//cpNode::cpNode(const cpNode& copy) : pnt(copy.pnt), index(copy.index),TEnode(copy.TEnode)
//{
//    
//}

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

edge* cpNode::getTE(edge* exception)
{
    for (int i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE() && edges[i] != exception)
        {
            if (edges[i]->getN2() == this)
            {
                edges[i]->flipDir();
            }
            return edges[i];
        }
    }
    return nullptr;
}

void cpNode::setTE() {TEnode = true;}

void cpNode::setIndex(int i) {index = i;}

edge* cpNode::getOtherTrailEdge(edge* current)
{
    for (int i=0; i<edges.size(); i++)
    {
        if (edges[i]->isTE() && edges[i] != current)
        {
            return edges[i];
        }
    }
    return nullptr;
}


