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
    bool upper; // Sheds wake panel from lower edge
    bool lower; // Sheds wake panel from upper edge
    bool lsFlag; // Lifting surface flag
    Eigen::Vector3d velocity;
    double Cp;
    
    int index; // Index in panel vector contained in geometry class.  Used for interpolating strength for wake panel influences.
    
    double srcSidePhi(const double &PN,const double &Al, const double &phiV,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s);
    Eigen::Vector3d srcSideV(const double &PN,const double &Al,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s,const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n);
    inline double pntSrcPhi(const double &PJK);
    inline Eigen::Vector3d pntSrcV(const Eigen::Vector3d &pjk);
    void tipVelocity();
    bool clusterTest(bodyPanel* other, double angle, const std::vector<bodyPanel*> &cluster);
    bool wingTipTest(bodyPanel* p);
    bool nearTrailingEdge();
//    bool neighborTests()
    std::vector<bodyPanel*> getBodyNeighbors();
    std::vector<bodyPanel*> gatherNeighbors(int nPanels);
    
public:
    bodyPanel(const Eigen::VectorXi &panelVertices,Eigen::MatrixXd* nodes, Eigen::Vector3d norm, int surfID,bool lsflag) : panel(panelVertices,nodes,norm,surfID), upper(false), lower(false), lsFlag(lsflag) {}
    
    
    bodyPanel(const bodyPanel &copy) : panel(copy), sourceStrength(copy.sourceStrength) {}
    
    void setNeighbors(panelOctree *oct, short normalMax);
    void setUpper() {upper = true;}
    void setLower() {lower = true;}
    void setIndex(int i) {index = i;}
    
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);
    
    void panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub);
    void panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub);
    
    void computeVelocity();
    void computeCp(double Vinf,double PG);
    Eigen::Vector3d computeMoments(const Eigen::Vector3d &cg);
    
    void setSigma(Eigen::Vector3d Vinf, double Vnorm)
    {
        sourceStrength = (-Vinf.dot(normal)+Vnorm);
    }
    
    void setMu(double dubStrength)
    {
        doubletStrength = dubStrength;
    }
    
    double getSigma() {return sourceStrength;}
    double getMu() {return doubletStrength;}
    bool isUpper() {return upper;}
    bool isLower() {return lower;}
    bool isLiftSurf() {return lsFlag;}
    int getIndex() {return index;}
    Eigen::Vector3d getGlobalV() {return velocity;}
    double getCp() {return Cp;}
};

#endif /* defined(__CPanel__bodyPanel__) */
