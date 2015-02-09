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
        n1->setTE();
        n2->setTE();
    }
    
    else if (bodyPans.size() == 2)
    {
        // Check for sharp edge without wake shed (i.e. vertical tail).  Used to start streamline tracing
        double angle = acos(bodyPans[0]->getNormal().dot(bodyPans[1]->getNormal()));
        if (angle > 5*M_PI/6 && bodyPans[0]->getID() == bodyPans[1]->getID())
        {
            TE = true;
            n1->setTE();
            n2->setTE();
            bodyPans[0]->setTEpanel(this);
        }
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

cpNode* edge::getOtherNode(cpNode* current)
{
    if (current == n1)
    {
        return n2;
    }
    else if (current == n2)
    {
        return n1;
    }
    else
    {
        return nullptr;
    }
}

double edge::length() {return (n2->getPnt()-n1->getPnt()).norm();}

std::vector<cpNode*> edge::getNodes()
{
    std::vector<cpNode*> ns;
    ns.push_back(n1);
    ns.push_back(n2);
    return ns;
}

Eigen::Vector3d edge::getVector()
{
    return (n2->getPnt()-n1->getPnt());
}

Eigen::Vector3d edge::getMidPoint()
{
    return (n1->getPnt()+0.5*getVector());
}

void edge::setNeighbors()
{
    if (bodyPans.size() == 2)
    {
        bodyPans[0]->addNeighbor(bodyPans[1]);
        bodyPans[1]->addNeighbor(bodyPans[0]);
    }
}

double edge::distToEdge(const Eigen::Vector3d &pnt)
{
    Eigen::Vector3d edgeVec = n2->getPnt()-n1->getPnt();
    Eigen::Vector3d pntVec = pnt-n1->getPnt();
    double dist = (pntVec.cross(edgeVec)).norm()/edgeVec.norm();
    return dist;
}