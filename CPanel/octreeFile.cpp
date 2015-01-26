//
//  octreeFile.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 11/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "octreeFile.h"

void octreeFile::writeFile(std::string filename,panelOctree* oct)
{
    std::vector<Eigen::Vector3d> centers;
    std::vector<Eigen::Vector3d> halfDs;
    std::vector<node<panel>*> nodes;
    nodes = oct->getNodes();
    
    for (int i=0; i<nodes.size(); i++)
    {
        centers.push_back(nodes[i]->getOrigin());
        halfDs.push_back(nodes[i]->getHalfDimension());
    }
    assert(centers.size() == halfDs.size());
    std::ofstream fid;
    fid.open(filename);
    if (fid.is_open())
    {
        fid << (nodes.size()+1) << "\n";
        for (int i=0; i<centers.size(); i++)
        {
            fid << centers[i](0) << "\t" << centers[i](1) << "\t" << centers[i](2) << "\n";
        }
        for (int i=0; i<halfDs.size(); i++)
        {
            fid << halfDs[i](0) << "\t" << halfDs[i](1) << "\t" << halfDs[i](2) << "\n";
        }
        fid.close();
    }
    else
    {
        std::cout << "Octree file could not be opened" << std::endl;
    }
}