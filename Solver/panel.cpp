//
//  panel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panel.h"

typedef Eigen::Vector3d         vector;
typedef Eigen::Vector3i         vertices;
typedef Eigen::Matrix3d         coordSys;


void panel::setGeom()
{
    longSide = 0;
    
    for (int i=0; i<verts.rows(); i++)
    {
        Eigen::Vector3d vec;
        double l;
        if (i != verts.rows()-1)
        {
            vec = nodes->row(verts(i+1))-nodes->row(verts(i));
        }
        else
        {
            vec = nodes->row(verts(0))-nodes->row(verts(i));
        }
        l = vec.norm();
        if (l>longSide)
        {
            longSide = l;
        }
    }
    
    if (verts.size() == 3)
    {
        Eigen::Vector3d p0,p1,p2;
        Eigen::Vector3d a,b;
        p0 = nodes->row(verts(0));
        p1 = nodes->row(verts(1));
        p2 = nodes->row(verts(2));
        a = p1-p0;
        b = p2-p0;
        center = (p0+p1+p2)/3;
        
        double theta = acos(a.dot(b)/(a.norm()*b.norm()));
        area = 0.5*a.norm()*b.norm()*sin(theta);
        
        normal = a.cross(b);
        normal.normalize();
    }
    else if (verts.size() == 4)
    {
        Eigen::Vector3d p0,p1,p2,p3,m1,m2;
        Eigen::Vector3d a,b,c,d,p,q;
        p0 = nodes->row(verts(0));
        p1 = nodes->row(verts(1));
        p2 = nodes->row(verts(2));
        p3 = nodes->row(verts(3));
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
    }
}

void panel::setNeighbors(panelOctree *oct)
{
    node<panel>* currentNode = oct->findNodeContainingMember(this);
    
    scanForNeighbors(currentNode,NULL);
    short maxNeighbors;
    if (TEpanel)
    {
        maxNeighbors = verts.size()+1;
        // Maximum number of neighbors on TE panel is the number of verts plus two.  This panel exists in the wing fuselage joint, with two neighbors on wing, two neighbors on fuselage, and either 1 (for tris) or 2 (for quads) in the wake.
    }
    else
    {
        maxNeighbors = verts.size();
    }
    
    while (currentNode != oct->getRootNode() && neighbors.size() < maxNeighbors)
    {
        scanForNeighbors(currentNode->getParent(),currentNode);
        if (TEpanel && ID > 10000 && neighbors.size()==verts.size()+1)
        {
            int wakePans = 0;
            for (int i=0; i<neighbors.size(); i++)
            {
                if (neighbors[i]->getID() > 10000)
                {
                    wakePans++;
                }
            }
            if (wakePans==verts.size()-2)
            {
                // Can only be wake panel in wing fuse joint so add one to max neighbors;
                maxNeighbors = verts.size()+2;
            }
        }
        currentNode = currentNode->getParent();
    }

}
    
void panel::scanForNeighbors(node<panel>* current, node<panel>* exception)
{
    std::vector<panel*> nodeMembers = current->getMembers(exception);
    for (int i=0; i<nodeMembers.size(); i++)
    {
        checkNeighbor(nodeMembers[i]);
    }
}

bool panel::isNeighbor(panel* other)
{
    Eigen::VectorXi otherVerts = other->getVerts();
    short count = 0;
    for (int i=0; i<verts.size(); i++)
    {
        for (int j=0; j<otherVerts.size(); j++)
        {
            if (nodes->row(verts(i)) == nodes->row(otherVerts(j)))
            {
                count++;
            }
            if (count == 2)
            {
                return true;
            }
        }
    }
    return false;
}

bool panel::neighborExists(panel* other)
{
    if (neighbors.size()>0)
    {
        for (int i=0; i<neighbors.size(); i++)
        {
            if (neighbors[i]==other)
            {
                return true;
            }
        }
    }

    return false;
}

void panel::checkNeighbor(panel* other)
{
    if (!neighborExists(other) && other != this) //Do not check and add panel if it is already a neighbor
    {
        if (isNeighbor(other))
        {
            addNeighbor(other);
            if (!other->isNeighbor(this))
            {
                other->addNeighbor(this);
            }
        }
    }
}

void panel::addNeighbor(panel* other)
{
    neighbors.push_back(other);
}

bool panel::isOnPanel(const Eigen::Vector3d &POI)
{
    Eigen::MatrixXd points(verts.rows()+1,3);
    for (int i=0; i<verts.rows(); i++)
    {
        points.row(i) = nodes->row(verts(i));
    }
    points.row(verts.rows()) = POI;
        
    convexHull hull(points,false);
    return (hull.getHull().size() == verts.size());
}

coordSys panel::getLocalSys()
{
    // Local Coordinate System
    // X : Points from center of panel to first vertex
    // Y : Normal crossed with X to obtain right hand coordinate system
    // Z : Normal to the panel
    coordSys local = Eigen::Matrix3d::Zero();
    local.row(0) = getUnitVector(center,nodes->row((verts(0))));
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

double panel::dubPhiInf(const Eigen::Vector3d &POI)
{
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm() < 0.00001)
    {
        return -0.5;
    }
    else if (std::abs(PN) < pow(10,-8))
    {
        return 0;
    }
    else if (pjk.norm()/longSide > 5)
    {
        return pntDubPhi(PN,pjk.norm());
    }
    else
    {
        double phi = 0;
        double Al;
        Eigen::Vector3d a,b,s;
        for (int i=0; i<verts.size(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            if (i!=verts.size()-1)
            {
                p1 = nodes->row(verts(i));
                p2 = nodes->row(verts(i+1));
            }
            else
            {
                p1 = nodes->row(verts(i));
                p2 = nodes->row(verts(0));
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
    if (pjk.norm() < 0.00001)
    {
        vel << 0,0,0;
        return vel;
    }
    else if (pjk.norm()/longSide > 5)
    {
        return pntDubV(local.row(2),pjk);
    }
    else
    {
        Eigen::Vector3d p1,p2,a,b,s;
        int i1,i2;
        for (int i=0; i<verts.size(); i++)
        {
            if (i!=verts.size()-1)
            {
                i1 = i;
                i2 = i+1;
            }
            else
            {
                i1 = i;
                i2 = 0;
            }
            p1 = nodes->row(verts(i1));
            p2 = nodes->row(verts(i2));
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
    double core = 0.1;
    return -(a.cross(b)*(a.norm()+b.norm()))/(a.norm()*b.norm()*((a.norm()*b.norm())+a.dot(b))+(pow(core*s.norm(),2)));
}

double panel::vortexPhi(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n)
{
    if (std::abs(PN) < pow(10,-10))
    {
        return 0;
    }
    else
    {
        double PA,PB,num,denom;
        
        PA = a.dot(l.cross(a.cross(s)));
        PB = PA-Al*s.dot(m);
        num = s.dot(m)*PN*(b.norm()*PA-a.norm()*PB);
        denom = PA*PB+pow(PN,2)*a.norm()*b.norm()*pow(s.dot(m),2);
        
        return atan2(num,denom);
    }
}
