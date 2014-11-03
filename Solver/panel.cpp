//
//  panel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panel.h"

typedef Eigen::Vector3d         vector;
typedef Eigen::Vector3d         point;
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
        point p0,p1,p2;
        vector a,b;
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
        point p0,p1,p2,p3,m1,m2;
        vector a,b,c,d,p,q;
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

void panel::setNeighbors(panelOctree *oct,bool wakePanel)
{
    node<panel>* currentNode = oct->findNodeContainingMember(this);
    
    scanForNeighbors(currentNode,NULL);
    short maxNeighbors;
    if (wakePanel)
    {
        maxNeighbors = verts.size()+2;
        // Maximum number of neighbors on wake panel is the number of verts plus two.  This panel exists in the wing fuselage joint, with two neighbors on wing, two neighbors on fuselage, and either 1 (for tris) or 2 (for quads) in the wake.
    }
    else
    {
        maxNeighbors = verts.size();
    }
    
    while (currentNode != oct->getRootNode() && neighbors.size() <= maxNeighbors)
    {
        scanForNeighbors(currentNode->getParent(),currentNode);
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

bool panel::isOnPanel(const point &POI)
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

void panel::checkNeighbor(panel* other)
{
    if (neighbors.size() <= verts.size() && !neighborExists(other) && other != this) //Do not check and add panel if it is already a neighbor or if maximum number of neighbors (4 for Tri, 5 for Quad) is already reached;
    {
        if (isNeighbor(other))
        {
            addNeighbor(other);
            other->addNeighbor(this);
        }
    }
}

void panel::addNeighbor(panel* other)
{
    neighbors.push_back(other);
}

coordSys panel::getLocalSys()
{
    // Local Coordinate System
    // X : Points from center of panel to first vertex
    // Y : Normal crossed with X to obtain right hand coordinate system
    // Z : Normal to the panel
    coordSys local;
    local.row(0) = getUnitVector(center,nodes->row((verts(0))));
    local.row(1) = normal.cross(local.row(0));
    local.row(2) = normal;
    
    return local;
}


vector panel::getUnitVector(const point &p1, const point &p2)
{
    //Returns unit vector pointing from p1 to p2;
    vector unit;
    for (int i=0; i<3; i++)
    {
        unit(i) = p2(i)-p1(i);
    }
    unit.normalize();
    return unit;
}

vector panel::transformCoordinates(const vector &toTransform, const coordSys &fromSys, const coordSys &toSys, const vector &translation)
{
    // Translation should be zeros if it is a vector transformation. It is need though for point transformations, such as the Point of Interest in the influence coefficient calculation.
    Eigen::Matrix3d transformCoeffs;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            transformCoeffs(i,j) = toSys.row(i).dot(fromSys.row(j));
        }
    }
    return transformCoeffs*(toTransform-translation);
}

double panel::doubletPhi(const double &mu, const point &POIglobal)
{
    
    // Establish panel coordinate system
    coordSys localSys = getLocalSys();
    coordSys globalSys;
    globalSys.setIdentity();
    
    // Transform Panel Vertices and Point of Interest to Local System
    vector POI = transformCoordinates(POIglobal,globalSys,localSys,center);
    
    if (POI.norm()/longSide > 5)
    {
        return area*POI(2)/(4*M_PI*pow(POI.norm(),3));
    }
    
    else
    {
        double phi=0;
        Eigen::MatrixXd vertsLocal(verts.rows(),3);
        for (int i=0; i<verts.rows(); i++)
        {
            vertsLocal.row(i) = transformCoordinates(nodes->row(verts(i)),globalSys,localSys,center);
        }
        if (d.rows() == 0)
        {
            setMD(vertsLocal);
        }
        
        Eigen::VectorXd r(verts.size()),e(verts.size()),h(verts.size());
        
        getREH(r,e,h,vertsLocal,POI);
        for (int i=0; i<vertsLocal.rows(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            double i1;
            double i2;
            if (i!=vertsLocal.rows()-1)
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(i+1);
                i1 = i;
                i2 = i+1;
            }
            else
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(0);
                i1 = i;
                i2 = 0;
            }
            phi = phi+1/(4*M_PI)*(atan2(m(i1)*e(i1)-h(i1),POI(2)*r(i1))-atan2(m(i1)*e(i2)-h(i2),POI(2)*r(i2)));
        }

        return phi;
    }
}

vector panel::doubletV(const double &mu, const point &POIglobal)
{
    
    // Establish panel coordinate system
    coordSys localSys = getLocalSys();
    coordSys globalSys;
    globalSys.setIdentity();
    
    // Transform Panel Vertices and Point of Interest to Local System
    vector POI = transformCoordinates(POIglobal,globalSys,localSys,center);
    
    vector vel;
    
    if (POI.norm()/longSide > 5)
    {
        vel(0) = 3*mu*area*POI(0)*POI(2)/(4*M_PI*pow(POI.norm(),5));
        vel(1) = 3*mu*area*POI(1)*POI(2)/(4*M_PI*pow(POI.norm(),5));
        vel(2) = -mu*area*(pow(POI(0),2)+pow(POI(1),2)-2*pow(POI(2),2))/(4*M_PI*pow(POI.norm(),5));
        return vel;
    }
    
    else
    {
        Eigen::MatrixXd vertsLocal(verts.rows(),3);
        for (int i=0; i<verts.rows(); i++)
        {
            vertsLocal.row(i) = transformCoordinates(nodes->row(verts(i)),globalSys,localSys,center);
        }
        if (d(0) == -10000)
        {
            setMD(vertsLocal);
        }
        Eigen::VectorXd r(verts.size()),e(verts.size()),h(verts.size());
        
        getREH(r,e,h,vertsLocal,POI);
        vector vTerms;
        for (int i=0; i<vertsLocal.rows(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            double i1;
            double i2;
            if (i!=vertsLocal.rows()-1)
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(i+1);
                i1 = i;
                i2 = i+1;
            }
            else
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(0);
                i1 = i;
                i2 = 0;
            }
            double denom = r(i1)*r(i2)*(r(i1)*r(i2)+((POI(0)-p1(0))*(POI(0)-p2(0))+(POI(1)-p1(1))*(POI(1)-p2(1))+pow(POI(2),2)));  // The + following r1*r2(r1*r2... is a - in Katz and Plotkin.  This does not yield velocities that are the same sign as the far field approximation.  Also, in the example programs in the back of Katz and Plotkin, it is a +.  Still waiting on getting a hold of the original Hess and Smith document from which it was drawn to confirm.
            
            vTerms(0) = vTerms(0)+(POI(2)*(p1(1)-p2(1))*(r(i1)+r(i2))/denom);
            vTerms(1) = vTerms(1)+(POI(2)*(p2(0)-p1(0))*(r(i1)+r(i2))/denom);
            vTerms(2) = vTerms(2)+(((POI(0)-p2(0))*(POI(1)-p1(1))-(POI(0)-p1(0))*(POI(1)-p2(1)))*(r(i1)+r(i2))/denom);
        }
        
        return -mu/(4*M_PI)*vTerms;;
    }
}

double panel::sourcePhi(const double &sigma, const point &POIglobal)
{
    // Establish panel coordinate system
    coordSys localSys = getLocalSys();
    coordSys globalSys;
    globalSys.setIdentity();
    
    // Transform Panel Vertices and Point of Interest to Local System
    vector POI = transformCoordinates(POIglobal,globalSys,localSys,center);
    
    if (POI.norm()/longSide > 5)
    {
        return -sigma*area/(4*M_PI*POI.norm());
    }
    
    else
    {
        Eigen::MatrixXd vertsLocal(verts.rows(),3);
        for (int i=0; i<verts.rows(); i++)
        {
            vertsLocal.row(i) = transformCoordinates(nodes->row(verts(i)),globalSys,localSys,center);
        }
        if (d(0) == -10000)
        {
            setMD(vertsLocal);
        }
        Eigen::VectorXd r(verts.size()),e(verts.size()),h(verts.size());
        
        getREH(r,e,h,vertsLocal,POI);
        double phiTerm1 = 0;
        double phiTerm2 = 0;
        for (int i=0; i<vertsLocal.rows(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            double i1;
            double i2;
            if (i!=vertsLocal.rows()-1)
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(i+1);
                i1 = i;
                i2 = i+1;
            }
            else
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(0);
                i1 = i;
                i2 = 0;
            }
            phiTerm1 = phiTerm1+((POI(0)-p1(0))*(p2(1)-p1(1))-(POI(1)-p1(1))*(p2(0)-p1(0)))/d(i1)*log((r(i1)+r(i2)+d(i1))/(r(i1)+r(i2)-d(i1)));
            
            phiTerm2 = phiTerm2+(atan2(m(i1)*e(i1)-h(i1),POI(2)*r(i1))-atan2(m(i1)*e(i2)-h(i2),POI(2)*r(i2)));
        }
        return sigma/(4*M_PI)*(phiTerm1-abs(POI(2))*phiTerm2);
        // Multiplied by negative one to account for traversing the perimeter of the element in a counter clockwise direction (per .tri format). Formulation from Hess and Smith is done based on a clockwise traverse of the perimeter.
    }
}

vector panel::sourceV(const double &sigma, const point &POIglobal)
{
    // Establish panel coordinate system
    coordSys localSys = getLocalSys();
    coordSys globalSys;
    globalSys.setIdentity();
    
    // Transform Panel Vertices and Point of Interest to Local System
    vector POI = transformCoordinates(POIglobal,globalSys,localSys,center);
    
    vector vel;
    
    if (POI.norm()/longSide > 5)
    {
        vel(0) = sigma*area*POI(0)/(4*M_PI*pow(POI.norm(),3));
        vel(1) = sigma*area*POI(1)/(4*M_PI*pow(POI.norm(),3));
        vel(2) = sigma*area*POI(2)/(4*M_PI*pow(POI.norm(),3));
        return vel;
    }
    
    else
    {
        Eigen::MatrixXd vertsLocal(verts.rows(),3);
        for (int i=0; i<verts.rows(); i++)
        {
            vertsLocal.row(i) = transformCoordinates(nodes->row(verts(i)),globalSys,localSys,center);
        }
        if (d(0) == -10000)
        {
            setMD(vertsLocal);
        }
        Eigen::VectorXd r(verts.size()),e(verts.size()),h(verts.size());
        
        getREH(r,e,h,vertsLocal,POI);
        vector vTerms;
        for (int i=0; i<vertsLocal.rows(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            double i1;
            double i2;
            if (i!=vertsLocal.rows()-1)
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(i+1);
                i1 = i;
                i2 = i+1;
            }
            else
            {
                p1 = vertsLocal.row(i);
                p2 = vertsLocal.row(0);
                i1 = i;
                i2 = 0;
            }
            
            vTerms(0) = vTerms(0)+(p2(1)-p1(1))/d(i1)*log((r(i1)+r(i2)-d(i1))/(r(i1)+r(i2)+d(i1)));
            vTerms(1) = vTerms(1)+(p1(0)-p2(0))/d(i1)*log((r(i1)+r(i2)-d(i1))/(r(i1)+r(i2)+d(i1)));
            vTerms(2) = vTerms(2)+(atan2(m(i1)*e(i1)-h(i1),POI(2)*r(i1))-atan2(m(i1)*e(i2)-h(i2),POI(2)*r(i2)));
        }
        return -sigma/(4*M_PI)*vTerms;
        // Multiplied by negative one to account for traversing the perimeter of the element in a counter clockwise direction (per .tri format). Formulation from Hess and Smith is done based on a clockwise traverse of the perimeter.
    }
}

void panel::getREH(Eigen::VectorXd &r, Eigen::VectorXd &e, Eigen::VectorXd &h,const Eigen::MatrixXd &verts, const point POI)
{
    for (int i=0; i<verts.rows(); i++)
    {
        Eigen::Vector3d p1;
        Eigen::Vector3d p2;
        if (i!=verts.rows()-1)
        {
            p1 = verts.row(i);
            p2 = verts.row(i+1);
            
        }
        else
        {
            p1 = verts.row(i);
            p2 = verts.row(0);
        }
        r(i) = sqrt(pow(POI(0)-p1(0),2)+pow(POI(1)-p1(1),2)+pow(POI(2),2));
        e(i) = pow(POI(0)-p1(0),2)+pow(POI(2),2);
        h(i) = (POI(0)-p1(0))*(POI(1)-p1(1));
    }

}

void panel::setMD(const Eigen::MatrixXd &verts)
{
    d.resize(verts.rows());
    m.resize(verts.rows());
    for (int i=0; i<verts.rows(); i++)
    {
        if (i == verts.rows()-1)
        {
            d(i) = sqrt(pow(verts(0,0)-verts(i,0),2)+pow(verts(0,1)-verts(i,1),2));
            m(i) = (verts(0,1)-verts(i,1))/(verts(0,0)-verts(i,0));
        }
        else
        {
            d(i) = sqrt(pow(verts(i+1,0)-verts(i,0),2)+pow(verts(i+1,1)-verts(i,1),2));
            m(i) = (verts(i+1,1)-verts(i,1))/(verts(i+1,0)-verts(i,0));
        }
    }
}
