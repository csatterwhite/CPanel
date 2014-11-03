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
    double mu;
    bodyPanel* upperPan;
    bodyPanel* lowerPan;
    bool TEpanel;
    std::vector<wakeLine*>* wakeLines;
    std::vector<bodyPanel*> interpPans;
    double interpCoeff;
    
    void addWakeLine();
    void boundingLines(wakeLine* &w1, wakeLine* &w2);
    
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
    wakePanel(const Eigen::VectorXi &panelVertices, Eigen::MatrixXd* nodes,int surfID,std::vector<wakeLine*>* wLines) : panel(panelVertices,nodes,surfID), wakeLines(wLines), TEpanel(false), interpCoeff(-10000) {};
    
    
    void setTEpanel(Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes);
    void findParentPanels(panelOctree* oct);
    double influenceCoeffPhi(const Eigen::Vector3d &POIglobal, std::vector<bodyPanel*> &interpPans,double &interpCoeff);
    double influenceCoeffV(const Eigen::Vector3d &POIglobal, std::vector<bodyPanel*> &interpPans,double &interpCoeff);
    double influencePhi(const Eigen::Vector3d &POIglobal);
    Eigen::Vector3d influenceV(const Eigen::Vector3d &POIglobal);

    bool isTEpanel() {return TEpanel;}
    void setMu();
    double getMu();
};

#endif /* defined(__CPanel__wakePanel__) */
