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

class wakePanel : public panel
{
    bodyPanel* upperPan;
    bodyPanel* lowerPan;
    std::vector<wakeLine*>* wakeLines;
    
    void addWakeLine();
    
    struct compareX
    {
        inline bool operator()(panel* p1, panel* p2)
        {
            return (p1->getCenter()(0) < p2->getCenter()(0));
        }
    };
    struct compareY
    {
        inline bool operator()(panel* p1, panel* p2)
        {
            return (p1->getCenter()(1) < p2->getCenter()(1));
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
    wakePanel(const Eigen::VectorXi &panelVertices, Eigen::MatrixXd* nodes,int surfID,std::vector<wakeLine*>* wLines) : panel(panelVertices,nodes,surfID), wakeLines(wLines), TEpanel(false) {};
    
    
    void setTEpanel(Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes);
    void findParentPanels(panelOctree* oct);
    void interpPanels(std::vector<bodyPanel*> &interpP, double &interpC);
    double panelPhi(const Eigen::Vector3d &POI);
    Eigen::Vector3d panelV(const Eigen::Vector3d &POI);

    bool isTEpanel() {return TEpanel;}
    void setMu();
    bodyPanel* getUpper() {return upperPan;}
    bodyPanel* getLower() {return lowerPan;}
};

#endif /* defined(__CPanel__wakePanel__) */
