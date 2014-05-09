//
//  geometry.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "geometry.h"
#include <fstream>
#include <array>
#include <set>


void geometry::readGeom(std::string geom_file)
{
    std::ifstream fid;
    fid.open(geom_file);
    if (fid.is_open())
    {
        fid >> nNodes >> nTris;
        Eigen::MatrixXi connectivity(nTris,3);
        Eigen::VectorXi surfID(nTris);
        std::vector<int> surfTypes;
        nodes.resize(nNodes,3);
        
        // Read XYZ Locations of Nodes
        for (int i=0; i<nNodes; i++)
        {
            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
        }
        
        // Temporarily Store Connectivity
        for (int i=0; i<nTris; i++)
        {
            fid >> connectivity(i,0) >> connectivity(i,1) >> connectivity(i,2);
        }
        
        // Scan Surface IDs and collect Unique IDs
        for (int i=0; i<nTris; i++)
        {
            fid >> surfID(i);
        }
        createSurfaces(connectivity,surfID);
        createOctree();
    }
}

void geometry::createSurfaces(Eigen::MatrixXi connectivity, Eigen::VectorXi surfID)
{
    surface* surf;
    for (int i=0; i<nTris; i++)
    {
        if (i==0 || surfID(i)!=surfID(i-1))
        {
            surf = new surface(surfID(i));
            surfaces.push_back(surf);
        }
        surfaces.back()->addPanel(connectivity.row(i),nodes);
    }
}

void geometry::createOctree()
{
    std::vector<panel*> data;
    for (int i=0; i<surfaces.size(); i++)
    {
        for (int j=0; j<surfaces[i]->getPanels().size(); j++)
        {
            data.push_back(surfaces[i]->getPanels()[j]);
        }
    }
    pOctree.addData(data);
}
