//
//  liftingSurf.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "liftingSurf.h"

liftingSurf::liftingSurf(int surfID,geometry* geom) : surface(surfID,geom)
{
    wakeSurf = new wake;
}

void liftingSurf::addPanel(bodyPanel* bPan)
{
    panels.push_back(bPan);
}
void liftingSurf::addPanel(wakePanel* wPan)
{
    wakeSurf->addPanel(wPan);
}
std::vector<panel*> liftingSurf::getAllPanels()
{
    std::vector<panel*> allPans;
    std::vector<wakePanel*> wakePans = wakeSurf->getPanels();
    allPans.insert(allPans.end(),panels.begin(),panels.end());
    allPans.insert(allPans.end(),wakePans.begin(),wakePans.end());
    return allPans;
}

