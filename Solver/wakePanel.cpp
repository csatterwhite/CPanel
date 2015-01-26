//
//  wakePanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakePanel.h"

wakePanel::wakePanel(const Eigen::VectorXi &panelVertices, Eigen::MatrixXd* nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, int surfID, wake* parentWake) : panel(panelVertices,nodes,bezNorm,surfID), pEdges(pEdges), TEpanel(false), parentWake(parentWake)
{
    for (int i=0; i<pEdges.size(); i++)
    {
        pEdges[i]->addWakePan(this);
    }
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

void wakePanel::setParentPanels()
{
    std::vector<edge*> TEedges;
    for (int i=0; i<pEdges.size(); i++)
    {
        if (pEdges[i]->isTE())
        {
            TEedges.push_back(pEdges[i]);
        }
    }
    std::vector<bodyPanel*> parentPans;
    
    // If TEedges.size() > 2, panel is at wing body joint
    for (int i=0; i<TEedges.size(); i++)
    {
        parentPans = TEedges[i]->getBodyPans();
        if (parentPans[0]->getID() == ID-10000)
        {
            // Test for edge on parent lifting surface if panel is at wing body joint
            break;
        }
    }
    if (parentPans[0]->getCenter()(2) > parentPans[1]->getCenter()(2))
    {
        upperPan = parentPans[0];
        lowerPan = parentPans[1];
        upperPan->setUpper();
        lowerPan->setLower();
    }
    else if (parentPans[0]->getCenter()(2) < parentPans[1]->getCenter()(2))
    {
        upperPan = parentPans[1];
        lowerPan = parentPans[0];
        upperPan->setUpper();
        lowerPan->setLower();
    }
    else
    {
        // Panels somehow are coplanar. Sort by normal direction
        if (parentPans[0]->getNormal()(2) > parentPans[1]->getNormal()(2))
        {
            upperPan = parentPans[0];
            lowerPan = parentPans[1];
            upperPan->setUpper();
            lowerPan->setLower();
        }
        else if (parentPans[0]->getNormal()(2) < parentPans[1]->getNormal()(2))
        {
            upperPan = parentPans[1];
            lowerPan = parentPans[0];
            upperPan->setUpper();
            lowerPan->setLower();
        }
    }
    
    wakeLine* wLine = new wakeLine(upperPan,lowerPan,normal);
    parentWake->addWakeLine(wLine);
}

//wakePanel* wakePanel::makeVortexSheet()
//{
//    Eigen::Vector3d p1,p2,p3,p4,deltaVec;
//    Eigen::VectorXi sheetVerts = Eigen::VectorXi::Zero(4);
//    double length = 100;
//    
//    p1(1) = -1000000;
//    
//    // Find points on trailing edge;
//    Eigen::Vector3i neighbVerts = upperPan->getVerts();
//    bool breakFlag = false;
//    for (int i=0; i<verts.size(); i++)
//    {
//        for (int j=0; j<neighbVerts.size(); j++)
//        {
//            if (nodes->row(verts(i)) == nodes->row(neighbVerts(j)))
//            {
//                Eigen::Vector3d pnt = nodes->row(verts(i));
//                if (pnt(1) > p1(1))
//                {
//                    sheetVerts(1) = sheetVerts(0);
//                    p2 = p1;
//                    
//                    sheetVerts(0) = verts(i);
//                    p1 = pnt;
//                    
//                }
//                else
//                {
//                    sheetVerts(1) = verts(i);
//                    p2 = pnt;
//                    breakFlag = true;
//                    break;
//                }
//            }
//        }
//        if (breakFlag)
//        {
//            break;
//        }
//    }
//    sheetVerts(2) = (int)nodes->rows();
//    sheetVerts(3) = sheetVerts(2) + 1;
//    
//    deltaVec = length*normal.cross((p2-p1).normalized());
//    p3 = p2+deltaVec;
//    p4 = p1+deltaVec;
//    
//    nodes->conservativeResize(nodes->rows()+2, 3);
//    nodes->row(sheetVerts(2)) = p3;
//    nodes->row(sheetVerts(3)) = p4;
//    
//    wakePanel* sheet = new wakePanel(sheetVerts,nodes,normal,ID,parentWake);
//    sheet->setTEpanel();
//    sheet->setUpper(upperPan);
//    sheet->setLower(lowerPan);
//    return sheet;
//    
//    
//}