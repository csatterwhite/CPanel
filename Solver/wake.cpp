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
    w = new wakePanel(panelVertices,nodes,surfID,TEnodes,&wakeLines);
    wpanels.push_back(w);
}

