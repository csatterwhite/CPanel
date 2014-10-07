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
#include "InfluenceTerms.h"

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
    
    point center;
    vector normal;
    vertices verts;
    double area;
    short surfID;
    bool TEpanel;
    std::vector<panel*> neighbors; //Share two vertices
    bool neighborExists(panel* other);
    bool isOnPanel(const point &POI, const Eigen::MatrixXd &nodes);
    vector getUnitVector(const point &p1, const point &p2);
    coordSys getLocalSys(const Eigen::MatrixXd &nodes);
    void pointSource(const double &sigma, const point &POI, double &phi, Eigen::Vector3d &vel);
    void pointDoublet(const double &mu, const point &POI, double &phi, Eigen::Vector3d &vel);
    void panelSource(const double &sigma, const point &POI, const Eigen::MatrixXd &vertsLocal, const influenceTerms &terms, double &phi, Eigen::Vector3d &vel);
    void panelDoublet(const double &mu, const point &POI, const Eigen::MatrixXd &vertsLocal, const influenceTerms &terms, double &phi, Eigen::Vector3d &vel);
    
public:
    panel(short surfaceID) : surfID(surfaceID) , TEpanel(false) {};
    void setGeom(const vertices &panelVertices, const Eigen::MatrixXd &nodes);
    void checkNeighbor(panel* other);
    void addNeighbor(panel* other);
    vector transformCoordinates(const vector &toTransform, const coordSys &fromSys, const coordSys &toSys);
    void sourceInfluence(const double &sigma, const point &POIglobal, const Eigen::MatrixXd &nodes, double &phi, Eigen::Vector3d &vel);
    void doubletInfluence(const double &mu, const point &POIglobal, const Eigen::MatrixXd &nodes, double &phi, Eigen::Vector3d &vel);
    
    
    vector getCenter() const {return center;}
    vector getNormal() const {return normal;}
    vertices getVerts() const {return verts;}
    std::vector<panel*> getNeighbors() const {return neighbors;}
    short getID() const {return surfID;}
    void setTE() {TEpanel = true;}
    bool isTEpanel() {return TEpanel;}
    
};

#endif /* defined(__CPanel__panel__) */
