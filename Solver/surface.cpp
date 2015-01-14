//
//  surface.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "surface.h"

//==== Destructor ====//
surface::~surface()
{
    for (int i=0; i<panels.size(); i++)
    {
        delete panels[i];
    }
}

void surface::addPanel(bodyPanel* bPan)
{
    panels.push_back(bPan);
}

void surface::setNeighbors(panelOctree* oct)
{
    for (int i=0; i<panels.size(); i++)
    {
        panels[i]->setNeighbors(oct, panels[i]->getVerts().size());
    }
}

