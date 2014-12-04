//
//  wake.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wake__
#define __CPanel__wake__

#include <stdio.h>
#include <vector>
#include "wakePanel.h"
#include "wakeLine.h"

class wake
{
    typedef Eigen::VectorXi vertices;
    
    std::vector<wakePanel*> wpanels;
    std::vector<wakeLine*> wakeLines;
    Eigen::MatrixXd* nodes;
    Eigen::Vector2d xlim;
    Eigen::Vector2d ylim;
    
    void setWakeDimensions();
    short edgePanel(wakePanel* p);
    
public:
    wake(Eigen::MatrixXd* nodes) : nodes(nodes) {}
    
    ~wake()
    {
        for (int i=0; i<wpanels.size(); i++)
        {
            delete wpanels[i];
        }
        for (int i=0; i<wakeLines.size(); i++)
        {
            delete wakeLines[i];
        }
    }
    
    wake(const wake& copy)
    {
        for (int i=0; i<copy.wpanels.size(); i++)
        {
            wpanels[i] = new wakePanel(*copy.wpanels[i]);
        }
    }
    
    void addPanel(const Eigen::VectorXi &panelVertices,Eigen::MatrixXd* nodes,Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes,int surfID);
    
    
    std::vector<wakePanel*> getPanels() const {return wpanels;}
    
    void setNeighbors(panelOctree* oct);
    
};

#endif /* defined(__CPanel__wake__) */
