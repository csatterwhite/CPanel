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
    for (int i=0; i<3; i++)
    {
        verts(i) = panelVertices(i);
    }
    setCenter(nodes);
    setNormal(nodes);
}

void panel::setCenter(const Eigen::MatrixXd &nodes)
{
    for (int i=0; i<3; i++)
    {
        double sum = 0;
        for (int j=0; j<3; j++)
        {
            sum += nodes(verts(j),i);
        }
        center(i) = sum/3;
    }
}

void panel::setNormal(const Eigen::MatrixXd &nodes)
{
    vector v01; //Unit vector from point 0 to 1
    vector v02; //Unit vector from point 0 to 2
    for (int i=0; i<3; i++)
    {
        v01(i) = nodes(verts(1),i)-nodes(verts(0),i);
        v02(i) = nodes(verts(2),i)-nodes(verts(0),i);
    }
    v01.normalize();
    v02.normalize();
    normal = v01.cross(v02);
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

bool panel::isOnPanel(std::vector<point> points, const point &POI)
{
    points.push_back(POI);
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

void panel::influenceTerms(const long &nVerts, const Eigen::MatrixXd &vertsLocal, const point &POI, const coordSys &localSys, const coordSys &globalSys, Eigen::VectorXd &d, Eigen::VectorXd &m, Eigen::VectorXd &r, Eigen::VectorXd &e, Eigen::VectorXd &h)
{
    // Calculates terms needed to compute source and doublet influence. Ref. Katz and Plotkin
    
    for (int i=0; i<nVerts; i++)
    {
        Eigen::Vector3d p1;
        Eigen::Vector3d p2;
        if (i!=verts.rows()-1)
        {
            p1 = vertsLocal.row(i);
            p2 = vertsLocal.row(i+1);
            
        }
        else
        {
            p1 = vertsLocal.row(i);
            p2 = vertsLocal.row(0);
        }
        d(i) = sqrt(pow(p2(0)-p1(0),2)+pow(p2(1)-p1(1),2));
        m(i) = (p2(1)-p1(1))/(p2(0)-p1(0));
        r(i) = sqrt(pow(POI(0)-p1(0),2)+pow(POI(1)-p1(1),2)+pow(POI(2),2));
        e(i) = pow(POI(0)-p1(0),2)+pow(POI(2),2);
        h(i) = (POI(0)-p1(0))*(POI(1)-p1(1));
    }
}

void panel::sourceInfluence(const double &sigma, const point &POIglobal, const Eigen::MatrixXd &nodes, vector &velocity, double &potential)
{
    // Establish panel coordinate system
    coordSys localSys = getLocalSys(nodes);
    coordSys globalSys;
    globalSys.setIdentity();
    // Transform Panel Vertices and Point of Interest to Local System
    vector POI = transformCoordinates(POIglobal,globalSys,localSys);
    
    long nVerts = verts.rows();
    Eigen::MatrixXd vertsLocal(nVerts,3);
    //  |x1 y1 z1| First Vertex
    //  |x2 y2 z2| Second Vertex
    //  |...     |
    //  |.       |
    //  |.       |
    //  |xN yN zN| Nth Vertex
    
    Eigen::VectorXd d(nVerts);
    Eigen::VectorXd m(nVerts);
    Eigen::VectorXd r(nVerts);
    Eigen::VectorXd e(nVerts);
    Eigen::VectorXd h(nVerts);
    
    for (int i=0; i<nVerts; i++)
    {
        vertsLocal.row(i) = transformCoordinates(nodes.row(verts(i)),globalSys,localSys);
    }
    
    influenceTerms(nVerts,vertsLocal,POI,localSys,globalSys,d,m,r,e,h);
    
    Eigen::MatrixXd velTerms(3,nVerts);
    
    for (int i=0; i<nVerts; i++)
    {
        Eigen::Vector3d p1;
        Eigen::Vector3d p2;
        double i1;
        double i2;
        if (i!=verts.rows()-1)
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
        velTerms(0,i) = (p2(1)-p1(1))/d(i1)*log((r(i1)+r(i2)-d(i1))/(r(i1)+r(i2)+d(i1)));
        velTerms(1,i) = (p1(0)-p2(0))/d(i1)*log((r(i1)+r(i2)-d(i1))/(r(i1)+r(i2)+d(i1)));
        velTerms(2,i) = atan2(m(i1)*e(i1)-h(i1),POI(2)*r(i1))-atan2(m(i1)*e(i2)-h(i2),POI(2)*r(i2));
    }
}
