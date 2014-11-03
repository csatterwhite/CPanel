//
//  wakePanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakePanel.h"

void wakePanel::findParentPanels(panelOctree* oct)
{
    setNeighbors(oct,true);
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
    
    std::sort(parentPanels.begin(),parentPanels.end(),compareX());
    std::sort(parentPanels.begin(),parentPanels.begin()+1,compareY());
    if (parentPanels.size() == 4)
    {
        // Special case where wake panel is at wing body joint and has two neighbors on wing and two neighbors on body. Filters neighbors on body so that wake panel parents are those on lifting surface.
        
        std::vector<bodyPanel*> temp;
        for (int i=0; i<parentPanels.size(); i++)
        {
            if (parentPanels[i]->getID() == ID-10000)
            {
                temp.push_back(parentPanels[i]);
            }
        }
        parentPanels = temp;
    }
    
    lowerPan = parentPanels[0];
    upperPan = parentPanels[1];
    lowerPan->setTEpanel();
    upperPan->setTEpanel();
    addWakeLine();
}

void wakePanel::setTEpanel(Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes)
{
    short count = 0;
    for (int i=0; i<verts.size(); i++)
    {
        if (TEnodes(verts(i)))
        {
            count++;
        }
        if (count == 2)
        {
            TEpanel = true;
            break;
        }
    }
}

void wakePanel::addWakeLine()
{
    wakeLine* wLine = new wakeLine(upperPan,lowerPan);
    wakeLines->push_back(wLine);
    std::sort(wakeLines->begin(),wakeLines->end(),compareLines());
}

double wakePanel::getMu()
{
    if (TEpanel)
    {
        mu = upperPan->getMu()-lowerPan->getMu();
    }
    return mu;
}

double wakePanel::influenceCoeffPhi(const Eigen::Vector3d &POIglobal, std::vector<bodyPanel*> &interpP, double &interpC)
{
    double IC = doubletPhi(1,POIglobal);
    if (interpCoeff == -10000)
    {
        wakeLine* wl1;
        wakeLine* wl2;
        boundingLines(wl1,wl2);
        interpCoeff = (center(1)-wl1->getY())/(wl2->getY()-wl1->getY());
        interpPans.push_back(wl1->getUpper());
        interpPans.push_back(wl1->getLower());
        interpPans.push_back(wl2->getUpper());
        interpPans.push_back(wl2->getLower());
    }
    interpC = interpCoeff;
    interpP = interpPans;
    return IC;
}

double wakePanel::influenceCoeffV(const Eigen::Vector3d &POIglobal, std::vector<bodyPanel*> &interpP, double &interpC)
{
    double IC = doubletV(1,POIglobal)(2);
    if (interpCoeff == -10000)
    {
        wakeLine* wl1;
        wakeLine* wl2;
        boundingLines(wl1,wl2);
        interpCoeff = (center(1)-wl1->getY())/(wl2->getY()-wl1->getY());
        interpPans.push_back(wl1->getUpper());
        interpPans.push_back(wl1->getLower());
        interpPans.push_back(wl2->getUpper());
        interpPans.push_back(wl2->getLower());
    }
    interpC = interpCoeff;
    interpP = interpPans;
    return IC;
}

double wakePanel::influencePhi(const Eigen::Vector3d &POIglobal)
{
    return doubletPhi(mu,POIglobal);
}

Eigen::Vector3d wakePanel::influenceV(const Eigen::Vector3d &POIglobal)
{
    return doubletV(mu,POIglobal);
}

void wakePanel::boundingLines(wakeLine* &w1, wakeLine* &w2)
{
    w1 = (*wakeLines)[0];
    w2 = (*wakeLines)[1];
    for (int i=0; i<wakeLines->size()-1; i++)
    {
        if (((*wakeLines)[i]->getY() <= center(1) && (*wakeLines)[i+1]->getY() > center(1)) || i == wakeLines->size()-2)
        {
            w1 = (*wakeLines)[i];
            w2 = (*wakeLines)[i+1];
        }
    }
}

void wakePanel::setMu()
{
    mu = (1-interpCoeff)*interpPans[0]->getMu() + (interpCoeff-1)*interpPans[1]->getMu() + interpCoeff*interpPans[2]->getMu() - interpCoeff*interpPans[3]->getMu();
}