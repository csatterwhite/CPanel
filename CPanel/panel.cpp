//
//  panel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panel.h"
#include "cpNode.h"
#include "edge.h"
#include "surface.h"

panel::panel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, int surfID) : ID(surfID), nodes(nodes), pEdges(pEdges), bezNormal(bezNorm)
{
    setGeom();
}

//panel::panel(const panel &copy) : ID(copy.ID), nodes(copy.nodes), pEdges(copy.pEdges)
//{
//    setGeom();
//}

void panel::setGeom()
{
    longSide = 0;
    
    for (int i=0; i<pEdges.size(); i++)
    {
        double l = pEdges[i]->length();
        if (l > longSide)
        {
            longSide = l;
        }
    }
    
    if (pEdges.size() == 3)
    {
        Eigen::Vector3d p0,p1,p2;
        Eigen::Vector3d a,b;
        p0 = nodes[0]->getPnt();
        p1 = nodes[1]->getPnt();
        p2 = nodes[2]->getPnt();
        a = p1-p0;
        b = p2-p0;
        center = (p0+p1+p2)/3;
        
        double theta = acos(a.dot(b)/(a.norm()*b.norm()));
        area = 0.5*a.norm()*b.norm()*sin(theta);
        normal = a.cross(b);
        normal.normalize();
        
        // If normals weren't included in input, set bezNormal to calculated normal
        if (bezNormal.isZero())
        {
            bezNormal = normal;
        }
        else
        {
            bezNormal.normalize();
        }
    }
    else if (pEdges.size() == 4)
    {
        Eigen::Vector3d p0,p1,p2,p3,m1,m2;
        Eigen::Vector3d a,b,c,d,p,q;
        p0 = nodes[0]->getPnt();
        p1 = nodes[1]->getPnt();
        p2 = nodes[2]->getPnt();
        p3 = nodes[3]->getPnt();
        a = p1-p0;
        b = p2-p1;
        c = p3-p2;
        d = p0-p3;
        p = b+c;
        q = a+b;
        m1 = p1+0.5*p;
        m2 = p0+0.5*q;
        center = 0.5*(m1+m2);
        area = 0.5*p.cross(q).norm();
        normal = a.cross(b);
        normal.normalize();
        
        // If normals weren't included in input, set bezNormal to calculated normal
        if (bezNormal.isZero())
        {
            bezNormal = normal;
        }
        else
        {
            bezNormal.normalize();
        }
    }
}

void panel::setPotential(Eigen::Vector3d Vinf)
{
    potential = Vinf.dot(center)-doubletStrength;
}

bool panel::inPanelProjection(const Eigen::Vector3d &POI, Eigen::Vector3d &projectedPnt)
{
    // Returns true if point is contained in extrusion of panel infinitely in normal direction
    Eigen::MatrixXd points(nodes.size()+1,3);
    std::vector<Eigen::Vector3d> nodesLocal;
    for (int i=0; i<nodes.size(); i++)
    {
        nodesLocal.push_back(global2local(nodes[i]->getPnt(), true));
        points.row(i) = nodesLocal[i];
    }
    points.row(nodes.size()) = global2local(POI,true);
    
    convexHull hull(points,true);
    
    if (hull.compareNodes(nodesLocal))
    {
        Eigen::Vector3d vec = POI-center;
        Eigen::Vector3d projVec = vec-(vec.dot(normal))*normal;
        projectedPnt = center + projVec;
//        projectedPnt = points.row(nodes.size());
//        projectedPnt(2) = 0; // Get point in panel plane.
//        projectedPnt = local2global(projectedPnt, true);
        return true;
    }
    
    projectedPnt = POI;
    return false;
}

Eigen::Matrix3d panel::getLocalSys()
{
    // Local Coordinate System
    // X : Points from center of panel to first vertex
    // Y : Normal crossed with X to obtain right hand coordinate system
    // Z : Normal to the panel
    Eigen::Matrix3d local = Eigen::Matrix3d::Zero();
    local.row(0) = getUnitVector(center,nodes[0]->getPnt());
    local.row(1) = normal.cross(local.row(0));
    local.row(2) = normal;
    
    return local;
}


Eigen::Vector3d panel::getUnitVector(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2)
{
    //Returns unit vector pointing from p1 to p2;
    Eigen::Vector3d unit;
    for (int i=0; i<3; i++)
    {
        unit(i) = p2(i)-p1(i);
    }
    unit.normalize();
    return unit;
}

Eigen::Vector3d panel::global2local(const Eigen::Vector3d &globalVec, bool translate)
{
    Eigen::Vector3d toTrans = globalVec;
    if (translate)
    {
        toTrans = globalVec-center;
    }
    return getLocalSys()*toTrans;
}

Eigen::Vector3d panel::local2global(const Eigen::Vector3d &localVec, bool translate)
{
    Eigen::Matrix3d transMat = getLocalSys();
    transMat.transposeInPlace();
    if (translate)
    {
        return transMat*localVec+center;
    }
    else
    {
        return transMat*localVec;
    }
    
}

double panel::dubPhiInf(const Eigen::Vector3d &POI)
{
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    
    if (pjk.norm() < 0.0000001)
    {
        return -0.5;
    }
    
    if (pjk.norm()/longSide > 5)
    {
        return pntDubPhi(PN,pjk.norm());
    }
    else
    {
        double phi = 0;
        double Al;
        Eigen::Vector3d a,b,s;
        for (int i=0; i<nodes.size(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            if (i!=nodes.size()-1)
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[i+1]->getPnt();
            }
            else
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[0]->getPnt();
            }
            a = POI-p1;
            b = POI-p2;
            s = p2-p1;
            Al = local.row(2).dot(s.cross(a));
            phi += vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
        }
        return phi/(4*M_PI);
    }
}


Eigen::Vector3d panel::dubVInf(const Eigen::Vector3d &POI)
{
    // VSAero doublet velocity influence formulation
    Eigen::Vector3d vel = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
//    if (pjk.norm() < 0.0000001)
//    {
//        vel << 0,0,0;
//        return vel;
//    }
    if (pjk.norm()/longSide > 5)
    {
        return pntDubV(local.row(2),pjk);
    }
    else
    {
        Eigen::Vector3d p1,p2,a,b,s;
        int i1,i2;
        for (int i=0; i<nodes.size(); i++)
        {
            if (i!=nodes.size()-1)
            {
                i1 = i;
                i2 = i+1;
            }
            else
            {
                i1 = i;
                i2 = 0;
            }
            p1 = nodes[i1]->getPnt();
            p2 = nodes[i2]->getPnt();
            a = POI-p1;
            b = POI-p2;
            s = p2-p1;
            
            vel += vortexV(a,b,s);
        }
        return vel/(4*M_PI);
    }
}

Eigen::Vector3d panel::vortexV(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s)
{
    double core = .05;
    return (a.cross(b)*(a.norm()+b.norm()))/(a.norm()*b.norm()*((a.norm()*b.norm())+a.dot(b))+(pow(core,2)));
}

double panel::vortexPhi(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n)
{
    double eps = pow(10, -15);
    double PA,PB,num,denom;
    
    PA = a.dot(l.cross(a.cross(s)));
    PB = PA-Al*s.dot(m);
    num = s.dot(m)*PN*(b.norm()*PA-a.norm()*PB);
    denom = PA*PB+pow(PN,2)*a.norm()*b.norm()*pow(s.dot(m),2);
    if (denom == 0 && std::abs(PN) < eps)
    {
        // Point is on edge.
        if (PN >= 0)
        {
            return 0.5*M_PI;
        }
        else
        {
            return -0.5*M_PI;
        }
    }
    return atan2(num,denom);
}

double panel::pntDubPhi(const double &PN, const double &PJK)
{
    return PN*area/(4*M_PI*pow(PJK,3));
}

Eigen::Vector3d panel::pntDubV(const Eigen::Vector3d n,const Eigen::Vector3d &pjk)
{
    return area*(3*pjk.dot(n)*pjk-pow(pjk.norm(),2)*n)/(4*M_PI*pow(pjk.norm(),5));
}

Eigen::VectorXi panel::getVerts()
{
    Eigen::VectorXi verts(nodes.size());
    for (int i=0; i<nodes.size(); i++)
    {
        verts(i) = nodes[i]->getIndex();
    }
    return verts;
}

std::vector<Eigen::Vector3d> panel::pntsAroundPnt(int nPnts,const Eigen::Vector3d &POI,double r)
{
    std::vector<Eigen::Vector3d> pnts;
    double theta;
    Eigen::Vector3d pnt;
    for (int i=0; i<nPnts; i++)
    {
        theta = (double)i/nPnts*(2*M_PI);
        pnt(0) = r*cos(theta);
        pnt(1) = r*sin(theta);
        pnt(2) = 0;
        pnt = local2global(pnt,true)+(POI-center);
        pnts.push_back(pnt);
    }
    
    return pnts;
}

Eigen::Vector3d panel::pntNearEdge(edge* e)
{
    Eigen::Vector3d pnt = center+0.95*(e->getMidPoint()-center);
    return pnt;
}

