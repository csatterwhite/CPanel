//
//  bodyPanel.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__bodyPanel__
#define __CPanel__bodyPanel__

#include <iostream>
#include "panel.h"

class bodyPanel : public panel
{
    double sourceStrength;
    bool TEpanel;
    
    double srcSidePhi(const double &PN,const double &Al, const double &phiV,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s);
    Eigen::Vector3d srcSideV(const double &PN,const double &Al,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s,const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n);
    inline double pntSrcPhi(const double &PJK);
    inline Eigen::Vector3d pntSrcV(const Eigen::Vector3d &pjk);
    
public:
    bodyPanel(const Eigen::VectorXi &panelVertices,Eigen::MatrixXd* nodes,int surfID) : panel(panelVertices,nodes,surfID), TEpanel(false) {}
    
    
    bodyPanel(const bodyPanel &copy) : panel(copy.verts,copy.nodes,copy.ID), sourceStrength(copy.sourceStrength), TEpanel(copy.TEpanel) {}
    
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);
    
    void panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub);
    void panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub);
    
    void setSigma(Eigen::Vector3d Vinf, double Vnorm)
    {
        sourceStrength = -Vinf.dot(normal)+Vnorm;
    }
    
    void setMu(double dubStrength)
    {
        doubletStrength = dubStrength;
    }
    
    void setTEpanel() {TEpanel = true;}
    
    double getSigma() {return sourceStrength;}
    double getMu() {return doubletStrength;}
    bool isTEpanel() {return TEpanel;}
};

#endif /* defined(__CPanel__bodyPanel__) */
