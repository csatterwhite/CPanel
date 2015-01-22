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

class wake;

class wakePanel : public panel
{
    bodyPanel* upperPan;
    bodyPanel* lowerPan;
    bool TEpanel;
    wake* parentWake;
    
    void neighborCases(const short &normalMax, short &absMax);
    
    struct compareX
    {
        inline bool operator()(panel* p1, panel* p2)
        {
            return (p1->getCenter()(0) < p2->getCenter()(0));
        }
    };
    struct compareZ
    {
        inline bool operator()(panel* p1, panel* p2)
        {
            return (p1->getCenter()(2) < p2->getCenter()(2));
        }
    };
    struct compareLines
    {
        inline bool operator()(wakeLine* w1, wakeLine* w2)
        {
            return (w1->getY() < w2->getY());
        }
    };
    
public:
    wakePanel(const Eigen::VectorXi &panelVertices, Eigen::MatrixXd* nodes, Eigen::Vector3d norm, int surfID, wake* parentWake) : panel(panelVertices,nodes,norm,surfID), TEpanel(false), parentWake(parentWake) {};
    
    wakePanel(const wakePanel &copy) : panel(copy),  upperPan(copy.upperPan), lowerPan(copy.lowerPan), TEpanel(copy.TEpanel), parentWake(copy.parentWake) {}
    
    void setNeighbors(panelOctree *oct, short normalMax);
    void setTEpanel() {TEpanel = true;}
    void setUpper(bodyPanel* up) {upperPan = up;}
    void setLower(bodyPanel* lp) {lowerPan = lp;}
    void setParentPanels();
    wakePanel* makeVortexSheet();
    void interpPanels(std::vector<bodyPanel*> &interpP, double &interpC);
    void addWakeLine(std::vector<wakeLine*> &wakeLines);
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);

    void setMu();
    void setStrength();
    bodyPanel* getUpper() {return upperPan;}
    bodyPanel* getLower() {return lowerPan;}
    bool isTEpanel() {return TEpanel;}
};

#endif /* defined(__CPanel__wakePanel__) */
