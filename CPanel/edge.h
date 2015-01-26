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
#include "bodyPanel.h"
#include "wakePanel.h"

class bodyPanel;
class wakePanel;

class edge
{
    std::vector<int> nodeIndices;
    std::vector<bodyPanel*> bodyPans;
    std::vector<wakePanel*> wakePans;
    bool TE; //Edge is at surface-wake junction
    
    void checkTE();
    
public:
    edge(int i1,int i2);
    
    void addBodyPan(bodyPanel* b);
    void addWakePan(wakePanel* w);
    bool sameEdge(int i1, int i2);
    
    bool isTE() {return TE;}
    std::vector<int> getIndices() {return nodeIndices;}
    std::vector<bodyPanel*> getBodyPans() {return bodyPans;}
    std::vector<wakePanel*> getWakePans() {return wakePans;}
};

#endif /* defined(__CPanel__edge__) */
