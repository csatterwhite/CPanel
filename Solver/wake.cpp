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

void wake::setNeighbors(panelOctree* oct)
{
    setWakeDimensions();
    for (int i=0; i<wpanels.size(); i++)
    {
        short nVerts = wpanels[i]->getVerts().size();
        if (edgePanel(wpanels[i]) == 2)
        {
            // Panel is on edge of wake and will have two neighbors for tris and three for quads
            wpanels[i]->setNeighbors(oct,nVerts-1);
        }
        else if (edgePanel(wpanels[i]) == 3)
        {
            // Panel is in downstream corner and will have 1 neighbor for tris and 2 for quads
            wpanels[i]->setNeighbors(oct,nVerts-2);
        }
        else
        {
            wpanels[i]->setNeighbors(oct,nVerts);
        }
        if (wpanels[i]->isTEpanel())
        {
            wpanels[i]->setParentPanels();
        }
    }
}

void wake::setWakeDimensions()
{
    Eigen::VectorXi verts;
    verts = wpanels[0]->getVerts();
    xlim(0) = nodes->row(verts(0))(0);
    xlim(1) = xlim(0);
    ylim(0) = nodes->row(verts(0))(1);
    ylim(1) = ylim(0);
    Eigen::Vector3d pnt;
    for (int i=0; i<wpanels.size(); i++)
    {
        verts = wpanels[i]->getVerts();
        for (int j=0; j<verts.size(); j++)
        {
            pnt = nodes->row(verts(j));
            if (pnt(0) < xlim(0))
            {
                xlim(0) = pnt(0);
            }
            else if (pnt(0) > xlim(1))
            {
                xlim(1) = pnt(0);
            }
            
            if (pnt(1) < ylim(0))
            {
                ylim(0) = pnt(1);
            }
            else if (pnt(1) > ylim(1))
            {
                ylim(1) = pnt(1);
            }
        }
    }
}

short wake::edgePanel(wakePanel* p)
{
    // Returns number of verts on edge of wake. 2 indicates panel edge of wake. 3 indicates panel in corner at downstream side of wake;
    short count = 0;
    Eigen::VectorXi verts = p->getVerts();
    Eigen::Vector3d pnt;
    for (int i=0; i<verts.size(); i++)
    {
        pnt = nodes->row(verts(i));
        if (pnt(0) == xlim(1) || pnt(1) == ylim(0) || pnt(1) == ylim(1))
        {
            count++;
        }
    }
    
    return count;
}