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
    
    if (parentPanels.size() == 4)
    {
        std::sort(parentPanels.begin(),parentPanels.end(),[](const panel* p1,const panel* p2) {return p1->getCenter()(0) < p2->getCenter()(0);}); // Ensures lifting surface panels are parent panels if the wake panel is at the wing body joint
        std::vector<bodyPanel*> npp; //In case where wake panel is at wing body joint, the parent panels will be selected as the two neighbors on the lifting surface, but the other two neighbors on the body need to be set as upper and lower to avoid taking the derivative of the potential across the discontinuous wake later on.
        npp.push_back(parentPanels[2]);
        npp.push_back(parentPanels[3]);
        std::sort(npp.begin(),npp.end(),[](const panel* p1,const panel* p2) {return p1->getCenter()(2) < p2->getCenter()(2);});
        npp[0]->setLower();
        npp[1]->setUpper();
    }
    
    std::vector<bodyPanel*> pp;
    pp.push_back(parentPanels[0]);
    pp.push_back(parentPanels[1]);
    std::sort(pp.begin(),pp.end(),[](const panel* p1,const panel* p2) {return p1->getCenter()(2) < p2->getCenter()(2);});
   
    lowerPan = pp[0];
    upperPan = pp[1];
    lowerPan->setLower();
    upperPan->setUpper();
}

void wakePanel::addWakeLine(std::vector<wakeLine*> &wakeLines)
{
    wakeLine* wLine = new wakeLine(upperPan,lowerPan,normal);
    wakeLines.push_back(wLine);
    std::sort(wakeLines.begin(),wakeLines.end(),compareLines());
}

void wakePanel::interpPanels(std::vector<bodyPanel*> &interpPans, double &interpCoeff)
{
    wakeLine* wl1 = nullptr;
    wakeLine* wl2 = nullptr;
    std::vector<wakeLine*> wakeLines = parentWake->getWakeLines();
    if (center(1) < wakeLines[1]->getY())
    {
        wl1 = wakeLines[0];
        wl2 = wakeLines[1];
    }
    else if (center(1) >= wakeLines.end()[-1]->getY())
    {
        wl1 = wakeLines.end()[-2]; //Second to last wakeline
        wl2 = wakeLines.end()[-1]; //Last wakeline
    }
    else
    {
        for (int i=1; i<wakeLines.size()-1; i++)
        {
            if ((wakeLines[i]->getY() <= center(1) && wakeLines[i+1]->getY() > center(1)))
            {
                wl1 = wakeLines[i];
                wl2 = wakeLines[i+1];
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

void wakePanel::setStrength()
{
    doubletStrength = upperPan->getMu()-lowerPan->getMu();
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
        }
        else
        {
            nls++;
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
            TEpanel = true;
        }
        else if (w < verts.size())
        {
            // Only two neighbors in wake means panel is shed from trailing edge and there will be one additional neighbor
            absMax = verts.size()+1;
            TEpanel = true;
        }
    }

}

horseshoeVortex* wakePanel::makeHorseshoe()
{
    Eigen::Vector3d p1,p2;

    p1(1) = 1000000;
    
    // Find points on trailing edge;
    Eigen::Vector3i neighbVerts = upperPan->getVerts();
    bool breakFlag = false;
    for (int i=0; i<verts.size(); i++)
    {
        for (int j=0; j<neighbVerts.size(); j++)
        {
            if (nodes->row(verts(i)) == nodes->row(neighbVerts(j)))
            {
                Eigen::Vector3d pnt = nodes->row(verts(i));
                if (pnt(1) < p1(1))
                {
                    p2 = p1;
                    p1 = pnt;
                    
                }
                else
                {
                    p2 = pnt;
                    breakFlag = true;
                    break;
                }
            }
        }
        if (breakFlag)
        {
            break;
        }
    }
    horseshoeVortex* h = new horseshoeVortex(p1,p2,upperPan,lowerPan);
    return h;
    
}

wakePanel* wakePanel::makeVortexSheet()
{
    Eigen::Vector3d p1,p2,p3,p4,deltaVec;
    Eigen::VectorXi sheetVerts = Eigen::VectorXi::Zero(4);
    double length = 100;
    
    p1(1) = -1000000;
    
    // Find points on trailing edge;
    Eigen::Vector3i neighbVerts = upperPan->getVerts();
    bool breakFlag = false;
    for (int i=0; i<verts.size(); i++)
    {
        for (int j=0; j<neighbVerts.size(); j++)
        {
            if (nodes->row(verts(i)) == nodes->row(neighbVerts(j)))
            {
                Eigen::Vector3d pnt = nodes->row(verts(i));
                if (pnt(1) > p1(1))
                {
                    sheetVerts(1) = sheetVerts(0);
                    p2 = p1;
                    
                    sheetVerts(0) = verts(i);
                    p1 = pnt;
                    
                }
                else
                {
                    sheetVerts(1) = verts(i);
                    p2 = pnt;
                    breakFlag = true;
                    break;
                }
            }
        }
        if (breakFlag)
        {
            break;
        }
    }
    sheetVerts(2) = (int)nodes->rows();
    sheetVerts(3) = sheetVerts(2) + 1;
    
    deltaVec = length*normal.cross((p2-p1).normalized());
    p3 = p2+deltaVec;
    p4 = p1+deltaVec;
    
    nodes->conservativeResize(nodes->rows()+2, 3);
    nodes->row(sheetVerts(2)) = p3;
    nodes->row(sheetVerts(3)) = p4;
    
    wakePanel* sheet = new wakePanel(sheetVerts,nodes,normal,ID,parentWake);
    sheet->setTEpanel();
    sheet->setUpper(upperPan);
    sheet->setLower(lowerPan);
    return sheet;
    
    
}