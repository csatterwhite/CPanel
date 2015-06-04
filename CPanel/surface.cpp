//
//  surface.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "surface.h"
#include "geometry.h"
#include "bodyStreamline.h"


surface::surface(const int &surfaceID,geometry* geom) : surfID(surfaceID), geom(geom), TEflag(false), LSflag(false) {}
//==== Destructor ====//
//surface::~surface()
//{
//    for (int i=0; i<panels.size(); i++)
//    {
//        delete panels[i];
//    }
//}

void surface::addPanel(bodyPanel* bPan)
{
    panels.push_back(bPan);
}

std::vector<edge*> surface::getTrailingEdges()
{
    if (!LSflag)
    {
        std::vector<edge*> empty;
        return empty;
    }
    else if (trailingEdges.size() > 0)
    {
        return trailingEdges;
    }
    
    edge* TE;
    
    for (int i=0; i<panels.size(); i++)
    {
        if (panels[i]->isTEpanel())
        {
            TE = panels[i]->getTrailingEdge();
            
            trailingEdges.push_back(TE);
        }
    }
    
    // Sort edges so that the streamlines can start tracing at one side of surface and progress, keeping track of neighboring streamlines.
    Eigen::Vector3d vec;
    vec = trailingEdges[0]->getVector();
    vec.normalize();
    
    if (vec(1) >= 0.01)
    {
        // Sort by Y direction (edges are primarily lateral)
        std::sort(trailingEdges.begin(), trailingEdges.end(), [](edge* e1,edge* e2) {return e1->getMidPoint()(1) < e2->getMidPoint()(1);});
    }
    else if (vec(1) < 0.01 && vec(2) >= 0.01)
    {
        // Sort by Z direction (edges are primarily vertical)
        std::sort(trailingEdges.begin(), trailingEdges.end(), [](edge* e1,edge* e2) {return e1->getMidPoint()(2) < e2->getMidPoint()(2);});
    }
    else
    {
        // Sort by X direction (edges are primarily streamwise)
        std::sort(trailingEdges.begin(), trailingEdges.end(), [](edge* e1,edge* e2) {return e1->getMidPoint()(2) < e2->getMidPoint()(2);});
    }
    return trailingEdges;
}

Eigen::Vector3d surface::rearStagnationPnt(const Eigen::Vector3d &Vinf, double PG, bodyPanel* &p)
{
    // Sort panels to get furthest downstream in front.
    std::sort(panels.begin(),panels.end(),[](bodyPanel* p1, bodyPanel* p2) {return (p1->getPotential() > p2->getPotential());});
    
    p = panels[0];
    bodyStreamline s(p->getCenter(),p,Vinf,PG,geom,20,true);
    
//    std::cout << s.getPnts().size() << std::endl;
    
//    std::vector<Eigen::Vector3d> pnts = s.getPnts();
//    for (int i=0; i<pnts.size(); i++)
//    {
//        std::cout << pnts[i](0) << "," << pnts[i](1) << "," << pnts[i](2) << ";" << std::endl;
//    }
    
    return s.getPnts().back();
}

std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> surface::getStreamlineStartPnts(const Eigen::Vector3d &Vinf, double PG)
{
    std::pair<Eigen::Vector3d,bodyPanel*> pair;
    bodyPanel* pStag;
    
    if (!LSflag)
    {
        std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> pntPairs;
        Eigen::Vector3d stagPnt = rearStagnationPnt(Vinf,PG,pStag);
        std::vector<Eigen::Vector3d> pnts = pStag->pntsAroundPnt(10,stagPnt,0.3*pStag->getLongSide());
        
        std::vector<bodyPanel*> closePanels = pStag->getRelatedPanels();
        Eigen::Vector3d projPnt;
//        std::sort(closePanels.begin(),closePanels.end(),[&](bodyPanel* p1,bodyPanel *p2) {return pStag->dist2Pan(p1) < pStag->dist2Pan(p2);});
        
//        std::ofstream fid;
//        fid.open("StreamlinePnts.txt");
        
        for (int i=0; i<pnts.size(); i++)
        {
            for (int j=0; j<closePanels.size(); j++)
            {
                if (closePanels[j]->inPanelProjection(pnts[i],projPnt))
                {
//                    Eigen::Vector3d v1,v2;
//                    v1 = closePanels[j]->getNormal();
//                    v2 = projPnt-closePanels[j]->getCenter();
//                    v1.normalize();
//                    v2.normalize();
//                    double dot = v1.dot(v2);
//                   
//                    std::cout << projPnt(0) << "\t" << projPnt(1) << "\t" << projPnt(2) << std::endl;
                    
                    pair = std::make_pair(projPnt,closePanels[j]);
                    pntPairs.push_back(pair);
                    break;
                }
            }
            
        }
//        fid.close();
        
        return pntPairs;
    }
    else
    {
        // Pairs 1 and 2 correspond to opposite sides of the sharp edge.
        std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> pairs1;
        std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> pairs2;


        std::vector<edge*> edges = getTrailingEdges();
        std::vector<bodyPanel*> pans;
        Eigen::Vector3d oldNorm,pnt;
        double theta1,theta2;
        for (int i=0; i<edges.size(); i++)
        {
            pans = edges[i]->getBodyPans();
            if (i == 0)
            {
                pair = std::make_pair(pans[0]->pntNearEdge(edges[i]), pans[0]);
                pairs1.push_back(pair);
                pair = std::make_pair(pans[1]->pntNearEdge(edges[i]),pans[1]);
                pairs2.push_back(pair);
                
                oldNorm = pans[0]->getNormal();
            }
            else
            {
                theta1 = pans[0]->getNormal().dot(oldNorm);
                theta2 = pans[1]->getNormal().dot(oldNorm);
                // Find which panel is on the same side as side 1. The dot product of that panels normal with the old side 1 normal would be larger than the dot product using the opposite sides normal.
                if (theta1 >= theta2)
                {
                    pair = std::make_pair(pans[0]->pntNearEdge(edges[i]), pans[0]);
                    pairs1.push_back(pair);
                    pair = std::make_pair(pans[1]->pntNearEdge(edges[i]),pans[1]);
                    pairs2.push_back(pair);
                    oldNorm = pans[0]->getNormal();
                }
                else
                {
                    pair = std::make_pair(pans[1]->pntNearEdge(edges[i]), pans[1]);
                    pairs1.push_back(pair);
                    pair = std::make_pair(pans[0]->pntNearEdge(edges[i]),pans[0]);
                    pairs2.push_back(pair);
                    oldNorm = pans[1]->getNormal();
                }
            }
        }
        pairs1.insert(pairs1.end(),pairs2.begin(),pairs2.end());
        return pairs1;
    }
}
                                                                                                 
                                                                                        
