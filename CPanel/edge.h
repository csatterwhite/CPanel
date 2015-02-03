//
//  edge.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/25/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__edge__
#define __CPanel__edge__

#include <stdio.h>
#include <Eigen/Dense>
#include <vector>
//#include "bodyPanel.h"
//#include "wakePanel.h"
//#include "cpNode.h"

class bodyPanel;
class wakePanel;
class cpNode;

class edge
{
    cpNode* n1;
    cpNode* n2;
    std::vector<bodyPanel*> bodyPans;
    std::vector<wakePanel*> wakePans;
    bool TE; //Edge is at surface-wake junction
    
    void checkTE();
    
public:
    edge(cpNode* n1,cpNode* n2);
    
    void addBodyPan(bodyPanel* b);
    void addWakePan(wakePanel* w);
    bool sameEdge(cpNode* node1, cpNode* node2);
    
    bool isTE();
    double length();
    std::vector<cpNode*> getNodes();
        
    std::vector<bodyPanel*> getBodyPans() {return bodyPans;}
    std::vector<wakePanel*> getWakePans() {return wakePans;}
    bodyPanel* getOtherBodyPan(bodyPanel* currentPan);
};

#endif /* defined(__CPanel__edge__) */
