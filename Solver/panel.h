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

class panel
{
    typedef Eigen::Vector3d         vector;
    typedef Eigen::Vector3i         vertices;
    
    vector center;
    vector normal;
    vertices verts;
    short surfID;
    bool TEpanel;
    std::vector<panel*> neighbors; //Share two vertices
    void setCenter(const Eigen::MatrixXd &nodes);
    void setNormal(const Eigen::MatrixXd &nodes);
    bool neighborExists(panel* other);
    
public:
    panel(short surfaceID) : surfID(surfaceID) , TEpanel(false) {};
    void setGeom(const Eigen::Vector3i &panelVertices, const Eigen::MatrixXd &nodes);
    void checkNeighbor(panel* other);
    void addNeighbor(panel* other);
    void setTE() {TEpanel = true;}
    
    vector getCenter() const {return center;}
    vector getNormal() const {return normal;}
    vertices getVerts() const {return verts;}
    std::vector<panel*> getNeighbors() const {return neighbors;}
    short getID() const {return surfID;}
    bool isTEpanel() {return TEpanel;}
    
};

#endif /* defined(__CPanel__panel__) */
