//
//  edge.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/25/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "edge.h"

edge::edge(int i1,int i2) : TE(false)
{
    nodeIndices.push_back(i1);
    nodeIndices.push_back(i2);
}

void edge::addBodyPan(bodyPanel* b)
{
    
    bodyPans.push_back(b);
    checkTE();
}

void edge::addWakePan(wakePanel* w)
{
    wakePans.push_back(w);
    checkTE();
}

void edge::checkTE()
{
    if (wakePans.size() == 1 && bodyPans.size() == 2)
    {
        TE = true;
        wakePans[0]->setTEpanel();
    }
}

bool edge::sameEdge(int i1, int i2)
{
    if (std::find(nodeIndices.begin(),nodeIndices.end(),i1) != nodeIndices.end() && std::find(nodeIndices.begin(),nodeIndices.end(),i2) != nodeIndices.end())
    {
        return true;
    }
    return false;
}