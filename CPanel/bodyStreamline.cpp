//
//  bodyStreamline.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/27/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "bodyStreamline.h"
#include "geometry.h"

bodyStreamline::bodyStreamline(Eigen::Vector3d startPnt,bodyPanel* startPan, const Eigen::Vector3d &Vinf,geometry* geom, int pntsPerPanel, bool marchFwd) :  pntsPerPanel(pntsPerPanel), geom(geom), Vinf(Vinf)
{
    Vmag = Vinf.norm();
    
    // Propogate streamline forwards or in reverse
    int marchDir;
    if (marchFwd)
    {
        marchDir = 1;
    }
    else
    {
        marchDir = -1;
    }
    
    // Allocate variables used in propogation
    Eigen::Vector3d pnt,vel,pntAbove,pntOnEdge,oldPnt,oldVel,normal,oldNorm;
    double tEdge,tNominal,pntPot;
    double eps = 0.0001;
    bool breakFlag = false;
    double maxAngle = 2*M_PI/3;
    edge* lastEdge = nullptr;
    
    // Set initial conditions
    pnt = startPnt;
    bodyPanel* currentPan = startPan;
    normal = currentPan->getNormal();
    tNominal = currentPan->getLongSide()/(pntsPerPanel*currentPan->getGlobalV().norm());
    std::vector<edge*> edges = currentPan->getEdges();
    oldVel = marchDir*Vinf;
    int tempPntsPerPanel = pntsPerPanel;

    // Euler integration to create streamline
    while (!breakFlag)
    {
        pntAbove = pnt + eps*normal;
        pntPot = geom->pntPotential(pntAbove,Vinf);
        vel = marchDir*currentPan->pntVelocity(pnt, pntPot);
        
        
        if (stagnationPnt(vel,oldVel,maxAngle))
        {
            break;
        }
        
        // Project velocity onto panel to correct for small error in CHTLS
        vel = vel-(vel.dot(normal))/normal.squaredNorm()*normal;
        
        oldPnt = pnt;
        oldVel = vel;
        pnts.push_back(pntAbove);
        velocities.push_back(vel);
        
        for (int i=0; i<edges.size(); i++)
        {
            if (edges[i] != lastEdge && edgeIntersection(edges[i], pnt, vel, tEdge, pntOnEdge))
            {
                if (tEdge > tNominal)
                {
                    // Next point is still on same panel
                    pnt = oldPnt+oldVel*tNominal;
                    pntAbove = pnt+eps*normal;
                    maxAngle = 0.75; // ~45 degrees
                }
                else
                {
                    // Next point is beyond panel edge so force point to be the point on the edge
                    oldNorm = currentPan->getNormal();
                    currentPan = edges[i]->getOtherBodyPan(currentPan);
                    normal = currentPan->getNormal();
                    pnt = pntOnEdge;
                    pntAbove = pnt+eps*(0.5*oldNorm+normal);
                    
                    lastEdge = edges[i];
                    edges =currentPan->getEdges();
                    
                    if (currentPan->getCp() > 0.95)
                    {
                        // Near stagnation point, decrease time step to avoid premature stagnation point detection
                        tempPntsPerPanel = 4*pntsPerPanel;
                    }
                    else
                    {
                        tempPntsPerPanel = pntsPerPanel;
                    }
                    
                    // For stagnation point detection, allow velocity vector to change ~30 degrees plus the change in panel normal
                    maxAngle = .75 + safeInvCos(oldNorm, currentPan->getNormal());
                    tNominal = currentPan->getLongSide()/(tempPntsPerPanel*currentPan->getGlobalV().norm());
                }
                break;
            }
        }
        if (pnt == oldPnt)
        {
            breakFlag = true;
        }
    }
}

bool bodyStreamline::edgeIntersection(edge* e,const Eigen::Vector3d &pnt, const Eigen::Vector3d &vel, double &dt, Eigen::Vector3d &pntOnEdge)
{
    //Algorithm for 3D line intersection can be found at mathworld.wolfram.com/Line-LineIntersection.html and comes from Goldman, R. "Intersection of Two Lines in Three-Space." Graphics Gems I. San Diego: Academic Press, p. 304, 1990.
    
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
        if (safeInvCos(pntOnEdge-p3,vel) < M_PI/2)
        {
            dt = (pntOnEdge-p3).norm()/vel.norm();
        }
        if (dt > 0)
        {
            return true;
        }
    }
    return false;
}

bool bodyStreamline::stagnationPnt(const Eigen::Vector3d vel, const Eigen::Vector3d &velOld, double maxAngle)
{
    double phi = safeInvCos(vel,velOld);
    if (phi > maxAngle || vel.norm() < 0.01*Vmag)
    {
        // Two velocity vectors are pointing at each other, indicating stagnation point has been reached.
        return true;
    }
    return false;
}

double bodyStreamline::safeInvCos(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2)
{
    double dot = v1.dot(v2)/(v1.norm()*v2.norm());
    if (dot > 1)
    {
        dot = 1;
    }
    else if (dot < -1)
    {
        dot = -1;
    }
    
    return acos(dot);

}
