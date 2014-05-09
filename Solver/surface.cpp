//
//  surface.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "surface.h"


void surface::addPanel(const vertices &verts, const coordinates &nodes)
{
    bodyPanel* b;
    wakePanel* w;
    
    if (surfID<=10000)
    {
        b = new bodyPanel(surfID);
        panels.push_back(b);
        panels.back()->setGeom(verts,nodes);
    }
    
    else
    {
        w = new wakePanel(surfID);
        panels.push_back(w);
        panels.back()->setGeom(verts,nodes);
    }
}