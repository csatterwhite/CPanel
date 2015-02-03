//
//  wakeLine.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakeLine.h"

wakeLine::wakeLine(bodyPanel* upper, bodyPanel* lower,Eigen::Vector3d normal) : upper(upper), lower(lower), normal(normal)
{
    setDimensions();
}

double wakeLine::getStrength()
{
    return upper->getMu()-lower->getMu();
}

void wakeLine::setDimensions()
{
    std::vector<edge*> es = upper->getEdges();
    std::vector<cpNode*> TEnodes;
    for (int i=0; i<es.size(); i++)
    {
        if (es[i]->isTE())
        {
            TEnodes = es[i]->getNodes();
        }
    }

    p1 = TEnodes[0]->getPnt();
    p2 = TEnodes[1]->getPnt();
    pMid = (p1+p2)/2;
}
