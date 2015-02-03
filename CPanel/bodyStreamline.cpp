//
//  bodyStreamline.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/27/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "bodyStreamline.h"

bodyStreamline::bodyStreamline(bodyPanel* startPan,std::vector<bodyPanel*>* bPanels,std::vector<wakePanel*>* wPanels, const Eigen::Vector3d &Vinf, int pntsPerPanel) : bPanels(bPanels), wPanels(wPanels), pntsPerPanel(pntsPerPanel)
{
    Eigen::Vector3d pnt,vel,pntAbove,pntOnEdge,oldPnt,oldVel,normal;
    oldVel = -Vinf;
    
    double tEdge; // time to reach next edge;
    double dt;
    int count = 1;
    double eps = 0.0000001;
    bool breakFlag = startPan->getStreamFlag();
    
    bodyPanel* currentPan = startPan;
    std::vector<edge*> edges = currentPan->getEdges();
    double pntPot;
    edge* lastEdge = nullptr;
    
    // Create start point
    if (startPan->isUpper() || startPan->isLower())
    {
        pnt = trailingEdgePnt(startPan);
//        pnt = startPan->getCenter();
    }
    else
    {
        pnt = startPan->getCenter();
    }

    while (!breakFlag)
    {
        pnt += 0.0001*currentPan->getNormal();

        normal = currentPan->getNormal();
        pntPot = pntPotential(pnt,Vinf,currentPan);
        vel = -currentPan->pntVelocity(pnt, pntPot); //Negative sign to propogate backwards
        
        // Project velocity onto panel to correct for small error in CHTLS
        vel = vel-(vel.dot(normal))/normal.squaredNorm()*normal;
        
        oldPnt = pnt;
        oldVel = vel;
        pntAbove = oldPnt + eps*normal;
        pnts.push_back(pntAbove);
        velocities.push_back(oldVel);
        potentials.push_back(pntPot);
        
        if (stagnationPnt(vel, oldVel))
        {
            break;
        }
        for (int i=0; i<edges.size(); i++)
        {
            if (edges[i] != lastEdge && edgeIntersection(edges[i], pnt, vel, tEdge, pntOnEdge))
            {
                if (count < pntsPerPanel-1)
                {
                    dt = count/(pntsPerPanel-1.0)*tEdge;
                    assert(dt != 0);
                    pnt = oldPnt+oldVel*dt;
                    count++;
                }
                else
                {
                    pnt = pntOnEdge;
                    assert(pnt != oldPnt);
                    currentPan->setStreamFlag();
                    currentPan = edges[i]->getOtherBodyPan(currentPan);
                    lastEdge = edges[i];
                    edges = currentPan->getEdges();
                    if (currentPan->getStreamFlag())
                    {
                        breakFlag = true;
                    }
                    count = 1;
                }
                break;
            }
        }
        if (pnt == oldPnt)
        {
            breakFlag = true;
        }
    }
//    for (int i=0; i<pnts.size(); i++)
//    {
//        std::cout << pnts[i](0) << "\t" << pnts[i](1) << "\t" << pnts[i](2) << std::endl;
//    }
//    std::cout << "\n" << std::endl;
//    for (int i=0; i<pnts.size(); i++)
//    {
//        std::cout << velocities[i](0) << "\t" << velocities[i](1) << "\t" << velocities[i](2) << std::endl;
//    }
//    std::cout << "\n" << std::endl;
//    for (int i=0; i<pnts.size(); i++)
//    {
//        std::cout << potentials[i] << std::endl;
//    }
}

double bodyStreamline::pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d Vinf, bodyPanel* currentPan)
{
    double pot = 0;
    for (int i=0; i<bPanels->size(); i++)
    {
        pot += (*bPanels)[i]->panelPhi(pnt);
    }
    for (int i=0; i<wPanels->size(); i++)
    {
        pot += (*wPanels)[i]->panelPhi(pnt);
    }
    pot += Vinf.dot(pnt);
    return pot;
}

Eigen::Vector3d bodyStreamline::trailingEdgePnt(bodyPanel* p)
{
    Eigen::Vector3d pnt,TEpnt;
    std::vector<cpNode*> TEnodes;
    for (int i=0; i<p->getEdges().size(); i++)
    {
        if (p->getEdges()[i]->isTE())
        {
            TEnodes = p->getEdges()[i]->getNodes();
            break;
        }
    }
    TEpnt = 0.5*(TEnodes[0]->getPnt()+TEnodes[1]->getPnt());
    // Move Point slightly upstream from trailing edge
    pnt = TEpnt + 0.1*(p->getCenter()-TEpnt);
    
    return pnt;
}

bool bodyStreamline::edgeIntersection(edge* e,const Eigen::Vector3d &pnt, const Eigen::Vector3d &vel, double &dt, Eigen::Vector3d &pntOnEdge)
{
    //Algorithm for 3D line intersection can be found at mathworld.wolfram.com/Line-LineIntersection.html and comes from Goldman, R. "Intersection of Two Lines in Three-Space." Graphics Gems I. Sand Diego: Academic Press, p. 304, 1990.
    
    // dt and pntOnEdge should only be used if true is returned, otherwise their behavior is undefined.
    
    Eigen::Vector3d p1,p2,p3,p4,a,c;
    double s;
    dt = -1;
    p1 = e->getNodes()[0]->getPnt();
    p2 = e->getNodes()[1]->getPnt();
    p3 = pnt;
    p4 = pnt + vel;
    
    a = p2-p1;
    // b in algorithm is equal to vel;
    c = p3-p1;
    s = c.cross(vel).dot(a.cross(vel))/a.cross(vel).squaredNorm();
    if (s >= 0 && s <= 1)
    {
        pntOnEdge = p1+a*s;
        for (int i=0; i<3; i++)
        {
            if (vel(i) != 0)
            {
                dt = (pntOnEdge(i)-p3(i))/vel(i);
                break;
            }
        }
        if (dt > 0)
        {
            return true;
        }
    }
    return false;
}

bool bodyStreamline::stagnationPnt(const Eigen::Vector3d vel, const Eigen::Vector3d &velOld)
{
    double dot = vel.dot(velOld)/(vel.norm()*velOld.norm());

    if (dot > 1)
    {
        dot = 1;
    }
    else if (dot < -1)
    {
        dot = -1;
    }
    double phi = acos(dot);
    if (phi > 2*M_PI/3)
    {
        // Two velocity vectors are pointing at each other, indicating stagnation point has been reached.
        return true;
    }
    return false;
}
