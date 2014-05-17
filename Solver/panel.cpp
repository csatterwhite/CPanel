//
//  panel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panel.h"


void panel::setGeom(const Eigen::Vector3i &panelVertices, const Eigen::MatrixXd &nodes)
{
    for (int i=0; i<3; i++)
    {
        verts(i) = panelVertices(i);
    }
    setCenter(nodes);
    setNormal(nodes);
}

void panel::setCenter(const Eigen::MatrixXd &nodes)
{
    for (int i=0; i<3; i++)
    {
        double sum = 0;
        for (int j=0; j<3; j++)
        {
            sum += nodes(verts(j),i);
        }
        center(i) = sum/3;
    }
}

void panel::setNormal(const Eigen::MatrixXd &nodes)
{
    vector v01; //Unit vector from point 0 to 1
    vector v02; //Unit vector from point 0 to 2
    for (int i=0; i<3; i++)
    {
        v01(i) = nodes(verts(1),i)-nodes(verts(0),i);
        v02(i) = nodes(verts(2),i)-nodes(verts(0),i);
    }
    v01.normalize();
    v02.normalize();
    normal = v01.cross(v02);
}

bool panel::neighborExists(panel* other)
{
    if (neighbors.size()>0)
    {
        for (int i=0; i<neighbors.size(); i++)
        {
            if (neighbors[i]==other)
            {
                return true;
            }
        }
    }

    return false;
}

void panel::checkNeighbor(panel* other)
{
    if (neighbors.size() <= 4 && !neighborExists(other) && other != this) //Do not check and add panel if it is already a neighbor or if maximum number of neighbors (4) is already reached;
    {
        vertices otherVerts = other->getVerts();
        short count = 0;
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                if (verts(i) == otherVerts(j))
                {
                    count++;
                }
            }
        }
        if (count==2)
        {
            addNeighbor(other);
            other->addNeighbor(this);
        }
    }
}

void panel::addNeighbor(panel* other)
{
    neighbors.push_back(other);
}


