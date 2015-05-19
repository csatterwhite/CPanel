//
//  bodyStreamline.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/27/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "bodyStreamline.h"
#include "geometry.h"

bodyStreamline::bodyStreamline(Eigen::Vector3d startPnt,bodyPanel* startPan, const Eigen::Vector3d &Vinf, double PG, geometry* geom, int pntsPerPanel, bool marchFwd) :  pntsPerPanel(pntsPerPanel), geom(geom), Vinf(Vinf), PG(PG)
{
    Vmag = Vinf.norm();
    
    // Propogate streamline forwards or in reverse
    if (marchFwd)
    {
        marchDir = 1;
    }
    else
    {
        marchDir = -1;
    }
    
    // Allocate variables used in propogation
    Eigen::Vector3d pnt,vel,pntAbove,pntOnEdge;
    double tEdge, dt, pntPot;
    double eps = 0.00001;
    bool stagPnt = false;
    double angleTol = M_PI/4;
    double maxAngle = angleTol;
    edge* lastEdge = nullptr;
    int pntsLeft = pntsPerPanel;
    edge* e = nullptr;
    
    std::vector<bodyPanel*> possiblePans;

    pnt = startPnt;
    possiblePans.push_back(startPan);
    
    pntAbove = pnt+eps*possiblePans[0]->getNormal(); // Ensure streamline is visible on exterior of panel for visualization
    
    pntPot = geom->pntPotential(pntAbove, Vinf);
    
    int i = 0;
    maxAngle = angleTol;
    
    while (i < possiblePans.size())
    {
        e = edgeIntersection(possiblePans[i], pnt, pntPot, vel, tEdge, pntOnEdge, maxAngle, lastEdge,stagPnt);
        if (stagPnt)
        {
            break;
        }
        if (e != nullptr)
        {
            pnts.push_back(pntAbove);
            velocities.push_back(vel);
            pntsLeft--;

            if (pntsLeft == 1)
            {
                pnt = pntOnEdge;
                possiblePans = e->getBodyPans();
                pntAbove = pnt + eps*e->getNormal();
                
                pntsLeft = pntsPerPanel; // Reset counter for next panel
                
                maxAngle = 2*angleTol+safeInvCos(possiblePans[0]->getNormal(), possiblePans[1]->getNormal());
                lastEdge = e;
            }
            else
            {
                dt = tEdge/(pntsLeft);
                pnt += vel*dt;
                pntAbove = pnt+eps*possiblePans[i]->getNormal();
                
                if (possiblePans.size() > 1)
                {
                    // Remove other possible panel from vector
                    for (int j=0; j<possiblePans.size(); j++)
                    {
                        if (possiblePans[j] != possiblePans[i])
                        {
                            possiblePans.erase(possiblePans.begin()+j);
                        }
                    }
                }
                
                maxAngle = angleTol;
            }
            pntPot = geom->pntPotential(pntAbove, Vinf);
            
            i = 0;
            continue;
        }
        else
        {
            i++;
        }
    }

}

edge* bodyStreamline::edgeIntersection(bodyPanel* pan,const Eigen::Vector3d &pnt, double pntPot, Eigen::Vector3d &vel, double &dt, Eigen::Vector3d &pntOnEdge, double maxAngle,edge* lastEdge, bool &stagFlag)
{
    //Algorithm for 3D line intersection can be found at mathworld.wolfram.com/Line-LineIntersection.html and comes from Goldman, R. "Intersection of Two Lines in Three-Space." Graphics Gems I. San Diego: Academic Press, p. 304, 1990.
    
    // dt and pntOnEdge should only be used if a non-null edge is returned, otherwise their behavior is undefined.
    edge* e = nullptr;
    
    Eigen::Vector3d p1,p2,p3,p4,a,c;
    double s;
    dt = -1;
    std::vector<edge*> edges = pan->getEdges();
    
    vel = marchDir*pan->pntVelocity(pnt, pntPot, PG);
    vel = vel-(vel.dot(pan->getNormal()))*pan->getNormal();
    
    // Check for Stagnation Point
    if (velocities.size() > 0)
    {
        stagFlag = stagnationPnt(vel, velocities.back(), maxAngle);
        if (stagFlag)
        {
            return nullptr;
        }
    }
    
    for (int i=0; i<edges.size(); i++)
    {
        if (edges[i] == lastEdge)
        {
            continue;
        }
        p1 = edges[i]->getNodes()[0]->getPnt();
        p2 = edges[i]->getNodes()[1]->getPnt();
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
                e = edges[i];
                break;
            }
        }
    }
    return e;
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
