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


void panel::setGeom(const vertices &panelVertices, const Eigen::MatrixXd &nodes)
{
    verts.resize(panelVertices.rows());
    for (int i=0; i<panelVertices.rows(); i++)
    {
        verts(i) = panelVertices(i);
    }
    
    if (verts.size() == 3)
    {
        point p0,p1,p2;
        vector a,b;
        p0 = nodes.row(verts(0));
        p1 = nodes.row(verts(1));
        p2 = nodes.row(verts(2));
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
        p0 = nodes.row(verts(0));
        p1 = nodes.row(verts(1));
        p2 = nodes.row(verts(2));
        p3 = nodes.row(verts(3));
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

bool panel::isOnPanel(const point &POI, const Eigen::MatrixXd &nodes)
{
    Eigen::MatrixXd points(verts.rows()+1,3);
    for (int i=0; i<verts.rows(); i++)
    {
        points.row(i) = nodes.row(verts(i));
    }
    points.row(verts.rows()) = POI;
        
    convexHull hull(points,false);
    return (hull.getHull().size() == verts.size());
}

void panel::checkNeighbor(panel* other)
{
    if (neighbors.size() <= 4 && !neighborExists(other) && other != this) //Do not check and add panel if it is already a neighbor or if maximum number of neighbors (4) is already reached;
    {
        vertices otherVerts = other->getVerts();
        short count = 0;
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                if (verts(i) == otherVerts(j))
                {
                    count++;
                }
            }
        }
        if (count==2)
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

coordSys panel::getLocalSys(const Eigen::MatrixXd &nodes)
{
    // Local Coordinate System
    // X : Points from center of panel to first vertex
    // Y : Normal crossed with X to obtain right hand coordinate system
    // Z : Normal to the panel
    coordSys local;
    local.row(0) = getUnitVector(center,nodes.row((verts(0))));
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

vector panel::transformCoordinates(const vector &toTransform, const coordSys &fromSys, const coordSys &toSys)
{
    Eigen::Matrix3d transformCoeffs;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            transformCoeffs(i,j) = toSys.row(i).dot(fromSys.row(j));
        }
    }
    return transformCoeffs*toTransform;
}

void panel::sourceInfluence(const double &sigma, const point &POIglobal, const Eigen::MatrixXd &nodes, double &phi, Eigen::Vector3d &vel)
{
    // Establish panel coordinate system
    coordSys localSys = getLocalSys(nodes);
    coordSys globalSys;
    globalSys.setIdentity();
    // Transform Panel Vertices and Point of Interest to Local System
    vector POI = transformCoordinates(POIglobal,globalSys,localSys);
    
    double eps = 0.0001;
    
    if (POI.norm()/area > 4)
    {
        pointSource(sigma,POI,phi,vel);
    }
    else if (POI.norm()/area<eps)
    {
        phi = 0;
        if (POI(2)>=0)
        {
            vel(2) = 0.5*sigma;
        }
        else
        {
            vel(2) = -0.5*sigma;
        }
        vel(0) = 0;
        vel(1) = 0;
    }
    else
    {
        Eigen::MatrixXd vertsLocal(verts.rows(),3);
        for (int i=0; i<verts.rows(); i++)
        {
            vertsLocal.row(i) = transformCoordinates(nodes.row(verts(i)), globalSys, localSys);
        }
        influenceTerms terms(vertsLocal,POI);
        panelSource(sigma,POI,vertsLocal,terms,phi,vel);
        
        // Handles special case where z->0 and z goes to zero when physically it should not be zero on the panel. Still need to handle case where point is near the boundary of the panel and the velocities become infinite.
        if (POI(2)/area<eps)
        {
            if (isOnPanel(POI,nodes))
            {
                if (POI(2)>=0)
                {
                    vel(2) = sigma/2;
                }
                else
                {
                    vel(2) = -sigma/2;
                }
            }
            else
            {
                vel(2) = 0;
            }
        }
    }
}

void panel::doubletInfluence(const double &mu, const point &POIglobal, const Eigen::MatrixXd &nodes, double &phi, Eigen::Vector3d &vel)
{
    // Establish panel coordinate system
    coordSys localSys = getLocalSys(nodes);
    coordSys globalSys;
    globalSys.setIdentity();
    
    // Transform Panel Vertices and Point of Interest to Local System
    vector POI = transformCoordinates(POIglobal,globalSys,localSys);
    
    double eps = 0.0001;
    
    if (POI.norm()/area > 4)
    {
        pointDoublet(mu,POI,phi,vel);
    }
    else if (POI.norm()/area<eps)
    {
        if (POI(2) >= 0)
        {
            phi = -mu/2;
        }
        else
        {
            phi = mu/2;
        }
        vel(0) = 0;
        vel(1) = 0;
        vel(2) = 0;
    }
    else
    {
        Eigen::MatrixXd vertsLocal(verts.rows(),3);
        for (int i=0; i<verts.rows(); i++)
        {
            vertsLocal.row(i) = transformCoordinates(nodes.row(verts(i)), globalSys, localSys);
        }
        influenceTerms terms(vertsLocal,POI);
        panelDoublet(mu,POI,vertsLocal,terms,phi,vel);
    }
}

void panel::pointSource(const double &sigma, const point &POI, double &phi, Eigen::Vector3d &vel)
{
    phi = -sigma*area/(4*M_PI*POI.norm());
    vel(0) = sigma*area*POI(0)/(4*M_PI*pow(POI.norm(),3));
    vel(1) = sigma*area*POI(1)/(4*M_PI*pow(POI.norm(),3));
    vel(2) = sigma*area*POI(2)/(4*M_PI*pow(POI.norm(),3));
}

void panel::pointDoublet(const double &mu, const point &POI, double &phi, Eigen::Vector3d &vel)
{
    phi = -mu*area*POI(2)/(4*M_PI*pow(POI.norm(),3));
    vel(0) = 3*mu*area*POI(0)*POI(2)/(4*M_PI*pow(POI.norm(),5));
    vel(1) = 3*mu*area*POI(1)*POI(2)/(4*M_PI*pow(POI.norm(),5));
    vel(2) = -mu*area*(pow(POI(0),2)+pow(POI(1),2)-2*pow(POI(2),2))/(4*M_PI*pow(POI.norm(),5));
}

void panel::panelSource(const double &sigma, const point &POI, const Eigen::MatrixXd &vertsLocal, const influenceTerms &terms, double &phi, Eigen::Vector3d &vel)
{
    Eigen::VectorXd d = terms.d;
    Eigen::VectorXd m = terms.m;
    Eigen::VectorXd r = terms.r;
    Eigen::VectorXd e = terms.e;
    Eigen::VectorXd h = terms.h;
    
    double phiTerm1 = 0;
    double phiTerm2 = 0;
    Eigen::Vector3d vTerms = Eigen::Vector3d::Zero();
    
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
        
        vTerms(0) = vTerms(0)+(p2(1)-p1(1))/d(i1)*log((r(i1)+r(i2)-d(i1))/(r(i1)+r(i2)+d(i1)));
        vTerms(1) = vTerms(1)+(p1(0)-p2(0))/d(i1)*log((r(i1)+r(i2)-d(i1))/(r(i1)+r(i2)+d(i1)));
    }
    vTerms(2) = phiTerm2;
    
    phi = sigma/(4*M_PI)*phiTerm1-abs(POI(2))*phiTerm2;
    vel = -sigma/(4*M_PI)*vTerms;
    // Multiplied by negative one to account for traversing the perimeter of the element in a counter clockwise direction (per .tri format). Formulation from Hess and Smith is done based on a clockwise traverse of the perimeter.
}

void panel::panelDoublet(const double &mu, const point &POI, const Eigen::MatrixXd &vertsLocal, const influenceTerms &terms, double &phi, Eigen::Vector3d &vel)
{
    Eigen::VectorXd d = terms.d;
    Eigen::VectorXd m = terms.m;
    Eigen::VectorXd r = terms.r;
    Eigen::VectorXd e = terms.e;
    Eigen::VectorXd h = terms.h;
    
    double phiTerm = 0;
    Eigen::Vector3d vTerms = Eigen::Vector3d::Zero();
    
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
        
        phiTerm = phiTerm+(atan2(m(i1)*e(i1)-h(i1),POI(2)*r(i1))-atan2(m(i1)*e(i2)-h(i2),POI(2)*r(i2)));
        
        double denom = r(i1)*r(i2)*(r(i1)*r(i2)+((POI(0)-p1(0))*(POI(0)-p2(0))+(POI(1)-p1(1))*(POI(1)-p2(1))+pow(POI(2),2)));  // The + following r1*r2(r1*r2... is a - in Katz and Plotkin.  This does not yield velocities that are the same sign as the far field approximation.  Also, in the example programs in the back of Katz and Plotkin, it is a +.  Still waiting on getting a hold of the original Hess and Smith document from which it was drawn to confirm.
        
        vTerms(0) = vTerms(0)+(POI(2)*(p1(1)-p2(1))*(r(i1)+r(i2))/denom);
        vTerms(1) = vTerms(1)+(POI(2)*(p2(0)-p1(0))*(r(i1)+r(i2))/denom);
        vTerms(2) = vTerms(2)+(((POI(0)-p2(0))*(POI(1)-p1(1))-(POI(0)-p1(0))*(POI(1)-p2(1)))*(r(i1)+r(i2))/denom);
    }
    
    phi = mu/(4*M_PI)*phiTerm;
    vel = -mu/(4*M_PI)*vTerms; // Multiplied by negative one to account for traversing the perimeter of the element in a counter clockwise direction (per .tri format). Formulation from Hess and Smith is done based on a clockwise traverse of the perimeter.
}


