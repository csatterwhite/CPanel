//
//  wake.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wake.h"


void wake::addPanel(const Eigen::VectorXi &panelVertices, Eigen::MatrixXd* nodes,Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes,int surfID)
{
    wakePanel* w;
    w = new wakePanel(panelVertices,nodes,surfID,&wakeLines);
    wpanels.push_back(w);
    wpanels.back()->setTEpanel(TEnodes);
}

void wake::setTEneighbors(panelOctree* oct)
{
    for (int i=0; i<wpanels.size(); i++)
    {
        if (wpanels[i]->isTEpanel())
        {
            wpanels[i]->findParentPanels(oct);
        }
    }
}
