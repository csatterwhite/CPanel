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
#include "math.h"
#include "convexHull.h"

class panel
{
    typedef Eigen::Vector3d         vector;
    typedef Eigen::Vector3d         point;
    typedef Eigen::Vector3i         vertices;
    typedef Eigen::Matrix3d         coordSys;
    
    //    coordSys
    //  |v1x v1y v1z| i.e. each axis is a row
    //  |v2x v2y v2z|
    //  |v3x v3y v3z|
    
    point center;
    vector normal;
    vertices verts;
    short surfID;
    bool TEpanel;
    std::vector<panel*> neighbors; //Share two vertices
    void setCenter(const Eigen::MatrixXd &nodes);
    void setNormal(const Eigen::MatrixXd &nodes);
    bool neighborExists(panel* other);
    bool isOnPanel(std::vector<point> points, const point &POI);
    vector getUnitVector(const point &p1, const point &p2);
    coordSys getLocalSys(const Eigen::MatrixXd &nodes);
    void influenceTerms(const long &nVerts, const Eigen::MatrixXd &nodes, const point &POIglobal, const coordSys &localSys, const coordSys &globalSys, Eigen::VectorXd &d, Eigen::VectorXd &m, Eigen::VectorXd &r, Eigen::VectorXd &e, Eigen::VectorXd &h);
    
public:
    panel(short surfaceID) : surfID(surfaceID) , TEpanel(false) {};
    void setGeom(const Eigen::Vector3i &panelVertices, const Eigen::MatrixXd &nodes);
    void checkNeighbor(panel* other);
    void addNeighbor(panel* other);
    vector transformCoordinates(const vector &toTransform, const coordSys &fromSys, const coordSys &toSys);
    void sourceInfluence(const double &sigma, const point &POIglobal, const Eigen::MatrixXd &nodes, vector &velocity, double &potential);
    
    
    vector getCenter() const {return center;}
    vector getNormal() const {return normal;}
    vertices getVerts() const {return verts;}
    std::vector<panel*> getNeighbors() const {return neighbors;}
    short getID() const {return surfID;}
    void setTE() {TEpanel = true;}
    bool isTEpanel() {return TEpanel;}
    
};

#endif /* defined(__CPanel__panel__) */
