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
#include <fstream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include "math.h"
#include "convexHull.h"
#include "panelOctree.h"
#include "chtlsnd.h"
//#include "cpNode.h"
//#include "edge.h"

class panelOctree;
class edge;
class cpNode;
class surface;

class panel
{    
protected:
    std::vector<cpNode*> nodes;
    std::vector<edge*> pEdges;
    Eigen::Vector3d center;
    Eigen::Vector3d normal;
    Eigen::Vector3d bezNormal; //Used in derivative calculation
    double area;
    double longSide;

    double doubletStrength;
    double potential;
    Eigen::Vector3d velocity;
    double Cp;
    int ID;
    
    Eigen::Vector3d global2local(const Eigen::Vector3d &globalVec,bool translate);
    Eigen::Vector3d local2global(const Eigen::Vector3d &localVec,bool translate);
    
    Eigen::Vector3d vortexV(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s);
    double vortexPhi(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n);
    Eigen::Vector3d getUnitVector(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);
    Eigen::Matrix3d getLocalSys();
    
    double pntDubPhi(const double &PN, const double &PJK);
    
    Eigen::Vector3d pntDubV(const Eigen::Vector3d n,const Eigen::Vector3d &pjk);
    
public:
    panel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm,int surfID);
    
    virtual ~panel() {}
    
//    panel(const panel &copy);
    
    void setGeom();
        
    void setPotential(Eigen::Vector3d Vinf);
    
    bool inPanelProjection(const Eigen::Vector3d &POI, Eigen::Vector3d &projectedPnt);
    
    double dubPhiInf(const Eigen::Vector3d &POI);
    Eigen::Vector3d dubVInf(const Eigen::Vector3d &POI);
    virtual double panelPhi(const Eigen::Vector3d &POI) = 0;
    virtual Eigen::Vector3d panelV(const Eigen::Vector3d &POI) = 0;
    
    std::vector<Eigen::Vector3d> pntsAroundPnt(int nPnts,const Eigen::Vector3d &POI, double r);
    Eigen::Vector3d pntNearEdge(edge* e);
    
    int getID() {return ID;}
    Eigen::Vector3d getCenter() const {return center;}
    Eigen::Vector3d getNormal() const {return normal;}
    Eigen::Vector3d getBezNormal() const {return bezNormal;}
    std::vector<cpNode*> getNodes() {return nodes;}
    Eigen::VectorXi getVerts();
    std::vector<edge*> getEdges() {return pEdges;}
    double getLongSide() {return longSide;}
    double getArea() {return area;}
    double getMu() {return doubletStrength;}
    double getPotential() {return potential;}
    
};

#endif /* defined(__CPanel__panel__) */
