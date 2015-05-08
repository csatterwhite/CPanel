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
#include "panelOctree.h"
#include "wakeLine.h"

class wake;
class bodyPanel;
class edge;

class wakePanel : public panel
{
    bodyPanel* upperPan;
    bodyPanel* lowerPan;
    bool TEpanel;
    wake* parentWake;
    
public:
    wakePanel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, wake* parentWake, int surfID);
    
    wakePanel(const wakePanel &copy);
    
    void setTEpanel();
    void setUpper(bodyPanel* up);
    void setLower(bodyPanel* lp);
    void setParentPanels(bodyPanel* upper, bodyPanel* lower);
    void setParentWake(wake* w) {parentWake = w;}
//    wakePanel* makeVortexSheet();
    void interpPanels(std::vector<bodyPanel*> &interpP, double &interpC);
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);

    void setMu();
    void setStrength();
    bodyPanel* getUpper() {return upperPan;}
    bodyPanel* getLower() {return lowerPan;}
    edge* getTE();
    bool isTEpanel() {return TEpanel;}
};

#endif /* defined(__CPanel__wakePanel__) */
