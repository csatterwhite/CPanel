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
//#include "edge.h"
//#include "cpNode.h"
//
//class edge;
class cpNode;

class bodyPanel : public panel
{
    surface* parentSurf;
    std::vector<bodyPanel*> neighbors;
    std::vector<bodyPanel*> cluster;
    int TSorder;
    double sourceStrength;
    bool upper; // Sheds wake panel from lower edge
    bool lower; // Sheds wake panel from upper edge
    bool TEpanel;
    edge* TE;
    bool tipFlag;
    bool streamFlag; // Surface Streamline crosses panel.
    Eigen::Vector3d velocity;
    double Cp;
    
    int index; // Index in panel vector contained in geometry class.  Used for interpolating strength for wake panel influences.
    
    double srcSidePhi(const double &PN,const double &Al, const double &phiV,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s);
    Eigen::Vector3d srcSideV(const double &PN,const double &Al,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s,const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n);
    inline double pntSrcPhi(const double &PJK);
    inline Eigen::Vector3d pntSrcV(const Eigen::Vector3d &pjk);
    bool clusterTest(bodyPanel* other, double angle,bool upFlag,bool lowFlag);
    bool wingTipTest();
    bool nearTrailingEdge();
    void setCluster(int dim, int bufferPanels);
    
public:
    bodyPanel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm,surface* parentSurf, int surfID);
    
//    bodyPanel(const bodyPanel &copy);
    
    void addNeighbor(bodyPanel* p);
    void setUpper();
    void setLower();
    void setTEpanel(edge* trailingEdge);
    void setIndex(int i);
    void setLSflag();
    void setTipFlag();
    
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);
    
    void panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub);
    void panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub);
    
    Eigen::Vector3d pntVelocity(const Eigen::Vector3d &pnt,double pntPot, double PG);
    
    double dist2Pan(bodyPanel* other);
    
    void computeVelocity(double PG);
    void computeCp(double Vinf);
    Eigen::Vector3d computeMoments(const Eigen::Vector3d &cg);
    
    void setSigma(Eigen::Vector3d Vinf, double Vnorm);
    
    void setMu(double dubStrength);
    
    void setStreamFlag();
    
    std::vector<bodyPanel*> getNeighbors() {return neighbors;}
    std::vector<bodyPanel*> getRelatedPanels();
    
    double getSigma() {return sourceStrength;}
    double getMu() {return doubletStrength;}
    bool isUpper() {return upper;}
    bool isLower() {return lower;}
    bool isTEpanel() {return TEpanel;}
    edge* getTrailingEdge() {return TE;}
    
    bool isTipPan() {return tipFlag;}
    bool getStreamFlag() {return streamFlag;}
    int getIndex() {return index;}
    Eigen::Vector3d getGlobalV() {return velocity;}
    double getCp() {return Cp;}
    
    std::vector<bodyPanel*> getCluster() {return cluster;}
};

#endif /* defined(__CPanel__bodyPanel__) */
