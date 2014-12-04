//
//  wakePanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakePanel.h"

void wakePanel::setParentPanels()
{
    std::vector<bodyPanel*> parentPanels;
    bodyPanel* ptr;
    for (int i=0; i<neighbors.size(); i++)
    {
        if (neighbors[i]->getID() != ID)
        {
            ptr = dynamic_cast<bodyPanel*>(neighbors[i]);
            parentPanels.push_back(ptr);
        }
    }
    
    std::sort(parentPanels.begin(),parentPanels.end(),compareX()); // Ensures lifting surface panels are parent panels if the wake panel is at the wing body joint
    std::sort(parentPanels.begin(),parentPanels.begin()+1,compareZ());
   
    lowerPan = parentPanels[0];
    upperPan = parentPanels[1];
    addWakeLine();
}

void wakePanel::addWakeLine()
{
    wakeLine* wLine = new wakeLine(upperPan,lowerPan);
    wakeLines->push_back(wLine);
    std::sort(wakeLines->begin(),wakeLines->end(),compareLines());
}

void wakePanel::interpPanels(std::vector<bodyPanel*> &interpPans, double &interpCoeff)
{
    wakeLine* wl1 = nullptr;
    wakeLine* wl2 = nullptr;
    if (center(1) < (*wakeLines)[1]->getY())
    {
        wl1 = (*wakeLines)[0];
        wl2 = (*wakeLines)[1];
    }
    else if (center(1) >= (*wakeLines).end()[-1]->getY())
    {
        wl1 = (*wakeLines).end()[-2]; //Second to last wakeline
        wl2 = (*wakeLines).end()[-1]; //Last wakeline
    }
    else
    {
        for (int i=1; i<wakeLines->size()-1; i++)
        {
            if (((*wakeLines)[i]->getY() <= center(1) && (*wakeLines)[i+1]->getY() > center(1)))
            {
                wl1 = (*wakeLines)[i];
                wl2 = (*wakeLines)[i+1];
            }
        }
    }
    assert(wl1 != nullptr || wl2 != nullptr);
    interpCoeff = (center(1)-wl1->getY())/(wl2->getY()-wl1->getY());
    
    interpPans[0] = wl1->getUpper();
    interpPans[1] = wl1->getLower();
    interpPans[2] = wl2->getUpper();
    interpPans[3] = wl2->getLower();
}


double wakePanel::panelPhi(const Eigen::Vector3d &POI)
{
    return doubletStrength*dubPhiInf(POI);
}

Eigen::Vector3d wakePanel::panelV(const Eigen::Vector3d &POI)
{
    return doubletStrength*dubVInf(POI);
}

void wakePanel::setMu()
{
    std::vector<bodyPanel*> interpPans(4);
    double interpCoeff;
    interpPanels(interpPans, interpCoeff);
    doubletStrength = (1-interpCoeff)*interpPans[0]->getMu() + (interpCoeff-1)*interpPans[1]->getMu() + interpCoeff*interpPans[2]->getMu() - interpCoeff*interpPans[3]->getMu();
}

void wakePanel::setNeighbors(panelOctree *oct, short normalMax)
{
    short absMax = normalMax;
    neighborCases(normalMax, absMax);
    node<panel>* currentNode = oct->findNodeContainingMember(this);
    node<panel>* exception = NULL;
    while (neighbors.size() < absMax)
    {
        scanForNeighbors(currentNode,exception);
        if (neighbors.size() >= absMax && normalMax > 1)
        {
            neighborCases(normalMax, absMax);
        }
        if (currentNode == oct->getRootNode())
        {
            break;
        }
        else
        {
            exception = currentNode;
            currentNode = currentNode->getParent();
        }
    }
}

void wakePanel::neighborCases(const short &normalMax, short &absMax)
{
    if (neighbors.size() < 2)
    {
        return;
    }
    short w = 0;
    short ls = 0;
    short nls = 0;
    for (int i=0; i<neighbors.size(); i++)
    {
        if (neighbors[i]->getID() == ID)
        {
            w++;
        }
        else if (neighbors[i]->getID() == ID-10000)
        {
            ls++;
            neighbors[i]->setTEpanel();
        }
        else
        {
            nls++;
            neighbors[i]->setTEpanel();
        }
    }
    if (ls > 0 || nls > 0)
    {
        TEpanel = true;
    }
    if (normalMax == verts.size()-1 && ls != 0)
    {
        // Wake panel is on outer edge of wake at the wing tip and will have three neighbors (2 body, 1 wake) for tris and four neighbors (2 body, 2 wake) for quads.
        absMax = verts.size();
    }
    if (normalMax == verts.size())
    {
        // Check special cases for wake panels shed off body
        if (w == verts.size())
        {
            // Only panels in the middle of the wake will have as many neighbors as vertices all in the wake
            return;
        }
        else if (ls > 0 && nls > 0)
        {
            // Corresponds to wake panel in the wing body joint
            absMax = verts.size()+2;
        }
        else if (w < verts.size())
        {
            // Only two neighbors in wake means panel is shed from trailing edge and there will be one additional neighbor
            absMax = verts.size()+1;
        }
        else
        {
            int dummy = 0;
        }
    }

}