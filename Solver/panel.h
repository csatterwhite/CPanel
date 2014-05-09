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
    std::vector<panel*> primaryNeighbors; //Share two vertices
    std::vector<panel*> secondaryNeighbors;  //Share one vertex
    void setCenter(const Eigen::MatrixXd &nodes);
    void setNormal(const Eigen::MatrixXd &nodes);
    // bool isNeighbor(panel* other);
    
public:
    panel(short surfaceID) : surfID(surfaceID) , TEpanel(false) {};
    void setGeom(const Eigen::Vector3i &panelVertices, const Eigen::MatrixXd &nodes);
    // void addNeighbor(panel* other);
    // void addPrimaryNeighbor(panel* other);
    // void addSecondaryNeighbor(panel* other);
    
    
    vector getCenter() const {return center;}
    vector getNormal() const {return normal;}
    vertices getVerts() const {return verts;}
    // std::vector<panel*> getPrimNeighbors() const {return primaryNeighbors;}
    // std::vector<panel*> getSecNeighbors() const {return secondaryNeighbors;}
    
};

#endif /* defined(__CPanel__panel__) */
