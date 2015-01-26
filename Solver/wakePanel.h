//
//  wakePanel.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wakePanel__
#define __CPanel__wakePanel__

#include <iostream>
#include "panel.h"
#include "bodyPanel.h"
#include "panelOctree.h"
#include "wakeLine.h"
#include "wake.h"
#include "edge.h"

class wake;
class edge;
class bodyPanel;

class wakePanel : public panel
{
    std::vector<edge*> pEdges;
    bodyPanel* upperPan;
    bodyPanel* lowerPan;
    bool TEpanel;
    wake* parentWake;
    
public:
    wakePanel(const Eigen::VectorXi &panelVertices, Eigen::MatrixXd* nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, int surfID, wake* parentWake);
    
    wakePanel(const wakePanel &copy) : panel(copy),  upperPan(copy.upperPan), lowerPan(copy.lowerPan), TEpanel(copy.TEpanel), parentWake(copy.parentWake) {}
    
    void setTEpanel() {TEpanel = true;}
    void setUpper(bodyPanel* up) {upperPan = up;}
    void setLower(bodyPanel* lp) {lowerPan = lp;}
    void setParentPanels();
//    wakePanel* makeVortexSheet();
    void interpPanels(std::vector<bodyPanel*> &interpP, double &interpC);
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);

    void setMu();
    void setStrength();
    bodyPanel* getUpper() {return upperPan;}
    bodyPanel* getLower() {return lowerPan;}
    bool isTEpanel() {return TEpanel;}
};

#endif /* defined(__CPanel__wakePanel__) */
