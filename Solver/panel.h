//
//  panel.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__panel__
#define __CPanel__panel__

#include <iostream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include "math.h"
#include "convexHull.h"
#include "panelOctree.h"

class panelOctree;

class panel
{
    typedef Eigen::Vector3d         vector;
    typedef Eigen::Vector3d         point;
    typedef Eigen::VectorXi         vertices;
    typedef Eigen::Matrix3d         coordSys;
    
    //    coordSys
    //  |v1x v1y v1z| i.e. each axis is a row
    //  |v2x v2y v2z|
    //  |v3x v3y v3z|
    
    double area;
    double longSide;
    Eigen::VectorXd d;
    Eigen::VectorXd m;

    bool neighborExists(panel* other);
    bool isNeighbor(panel* other);
    void scanForNeighbors(node<panel>* current, node<panel>* exception);
    bool isOnPanel(const point &POI);
    void checkNeighbor(panel* other);
    vector getUnitVector(const point &p1, const point &p2);
    coordSys getLocalSys();
    void getREH(Eigen::VectorXd &r, Eigen::VectorXd &e, Eigen::VectorXd &h,const Eigen::MatrixXd &verts, const point POI);
    void setMD(const Eigen::MatrixXd &verts);
    
protected:
    point center;
    vector normal;
    vertices verts;
    double potential;
    Eigen::MatrixXd* nodes;
    int ID;
    std::vector<panel*> neighbors; //Share two vertices
    
public:
    panel(const vertices &panelVertices,Eigen::MatrixXd* nodes,int surfID) : ID(surfID), verts(panelVertices), nodes(nodes)
    {
//        d.resize(verts.rows());
//        m.resize(verts.size());
//        d(0) = -10000; //Used as catch to set d and m if they haven't been set yet.
        setGeom();
        
    };
    
    virtual ~panel() {}
    
    panel(const panel &copy) : ID(copy.ID), verts(copy.verts), nodes(copy.nodes)
    {
        setGeom();
    }
    
    void setGeom();
    void setNeighbors(panelOctree *oct,bool wakePanel);
    void addNeighbor(panel* other);
    vector transformCoordinates(const vector &toTransform, const coordSys &fromSys, const coordSys &toSys, const vector &translation);
    double doubletPhi(const double &mu, const point &POIglobal);
    vector doubletV(const double &mu, const point &POIglobal);
    double sourcePhi(const double &sigma, const point &POIglobal);
    vector sourceV(const double &sigma, const point &POIglobal);
    
    int getID() {return ID;}
    vector getCenter() const {return center;}
    vector getNormal() const {return normal;}
    vertices getVerts() const {return verts;}
    std::vector<panel*> getNeighbors() const {return neighbors;}
    double getPotential() {return potential;}
};

#endif /* defined(__CPanel__panel__) */
