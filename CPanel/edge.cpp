//
//  edge.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/25/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "edge.h"
#include "bodyPanel.h"
#include "wakePanel.h"
#include "cpNode.h"

edge::edge(cpNode* n1,cpNode* n2) : n1(n1), n2(n2), TE(false)
{
    n1->addEdge(this);
    n2->addEdge(this);
}

void edge::addBodyPan(bodyPanel* b)
{
    
    bodyPans.push_back(b);
    checkTE();
}

void edge::addWakePan(wakePanel* w)
{
    wakePans.push_back(w);
    checkTE();
}

void edge::checkTE()
{
    if (wakePans.size() == 1 && bodyPans.size() == 2)
    {
        TE = true;
        wakePans[0]->setTEpanel();
    }
}

bool edge::isTE()  {return TE;}

bool edge::sameEdge(cpNode* node1, cpNode* node2)
{
    if ((node1 == n1 && node2 == n2) || (node1 == n2 && node2 == n1))
    {
        return true;
    }
    return false;
}

bodyPanel* edge::getOtherBodyPan(bodyPanel* currentPan)
{
    for (int i=0; i<2; i++)
    {
        if (bodyPans[i] != currentPan)
        {
            return bodyPans[i];
        }
    }
    return nullptr;
}

double edge::length() {return (n2->getPnt()-n1->getPnt()).norm();}

std::vector<cpNode*> edge::getNodes()
{
    std::vector<cpNode*> ns;
    ns.push_back(n1);
    ns.push_back(n2);
    return ns;
}