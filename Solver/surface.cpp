//
//  surface.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "surface.h"


surface::~surface()
{
    for (int i=0; i<panels.size(); i++)
    {
        delete panels[i];
    }
}

void surface::addPanel(const vertices &verts,Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes)
{
    bodyPanel* b;
    
    b = new bodyPanel(verts,nodes,surfID,TEnodes);
    panels.push_back(b);
}

void surface::setNeighbors(panelOctree* oct)
{
    for (int i=0; i<panels.size(); i++)
    {
        panels[i]->setNeighbors(oct, panels[i]->getVerts().size());
    }
}

