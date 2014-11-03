//
//  liftingSurf.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "liftingSurf.h"


void liftingSurf::addPanel(const Eigen::VectorXi &verts,Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes,int surfID)
{
    if (surfID <= 10000)
    {
        surface::addPanel(verts);
    }
    else
    {
        wakeSurf->addPanel(verts,nodes,TEnodes,surfID);
    }
}
std::vector<panel*> liftingSurf::getAllPanels()
{
    std::vector<panel*> allPans;
    std::vector<wakePanel*> wakePans = wakeSurf->getPanels();
    allPans.insert(allPans.end(),panels.begin(),panels.end());
    allPans.insert(allPans.end(),wakePans.begin(),wakePans.end());
    return allPans;
}

void liftingSurf::setTEneighbors(panelOctree* oct)
{
    wakeSurf->setTEneighbors(oct);
}