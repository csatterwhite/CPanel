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
protected:
    Eigen::Vector3d center;
    Eigen::Vector3d normal;
    Eigen::VectorXi verts;
    double area;
    double longSide;
    
    bool TEpanel;

    double doubletStrength;
    double potential;
    Eigen::MatrixXd* nodes;
    int ID;
    std::vector<panel*> neighbors;
    
    bool neighborExists(panel* other);
    bool isNeighbor(panel* other);
    void scanForNeighbors(node<panel>* current, node<panel>* exception);
    bool isOnPanel(const Eigen::Vector3d &POI);
    void checkNeighbor(panel* other);
    
    Eigen::Vector3d vortexV(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &s);
    double vortexPhi(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s, const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n);
    Eigen::Vector3d getUnitVector(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);
    Eigen::Matrix3d getLocalSys();
    
    inline double pntDubPhi(const double &PN, const double &PJK)
    {
        return PN*area/(4*M_PI*pow(PJK,3));
    }
    
    inline Eigen::Vector3d pntDubV(const Eigen::Vector3d n,const Eigen::Vector3d &pjk)
    {
        return -area*(3*pjk.dot(n)*pjk-pow(pjk.norm(),2)*n)/(4*M_PI*pow(pjk.norm(),5));
    }
    
public:
    panel(const Eigen::VectorXi &panelVertices,Eigen::MatrixXd* nodes,int surfID,Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes) : ID(surfID), verts(panelVertices), nodes(nodes), TEpanel(false)
    {
        setGeom();
    };
    
    virtual ~panel() {}
    
    panel(const panel &copy) : ID(copy.ID), verts(copy.verts), nodes(copy.nodes)
    {
        setGeom();
    }
    void setGeom();
    virtual void setNeighbors(panelOctree *oct, short normalMax) = 0;
    void setTEpanel()
    {
        TEpanel = true;
    }
    void setPotential(Eigen::Vector3d Vinf)
    {
        potential = Vinf.dot(center)-doubletStrength;
    }
    void addNeighbor(panel* other);
    Eigen::MatrixXd getLocalVerts();
    
    double dubPhiInf(const Eigen::Vector3d &POI);
    Eigen::Vector3d dubVInf(const Eigen::Vector3d &POI);
    virtual double panelPhi(const Eigen::Vector3d &POI) = 0;
    virtual Eigen::Vector3d panelV(const Eigen::Vector3d &POI) = 0;
    
    int getID() {return ID;}
    Eigen::Vector3d getCenter() const {return center;}
    Eigen::Vector3d getNormal() const {return normal;}
    Eigen::VectorXi getVerts() const {return verts;}
    double getArea() {return area;}
    bool isTEpanel() {return TEpanel;}
    std::vector<panel*> getNeighbors() const {return neighbors;}
    double getMu() {return doubletStrength;}
    double getPotential() {return potential;}
};

#endif /* defined(__CPanel__panel__) */
