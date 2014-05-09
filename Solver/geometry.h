//
//  geometry.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__geometry__
#define __CPanel__geometry__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <unordered_map>
#include "panelOctree.h"
#include "surface.h"
#include "wakePanel.h"
#include "bodyPanel.h"

class geometry
{
    std::vector<surface*> surfaces;
    panelOctree pOctree;
    Eigen::MatrixXd nodes;
    short nNodes;
    short nTris;

    void createSurfaces(Eigen::MatrixXi connectivity, Eigen::VectorXi surfID);
    void createOctree();
    // void setTEPanels();
    
public:
    geometry() {}
    
    ~geometry()
    {
        for (int i=0; i<surfaces.size(); i++)
        {
            delete surfaces[i];
        }
    }
    
    geometry(const geometry& copy) : pOctree(copy.pOctree), nodes(copy.nodes), nNodes(copy.nNodes), nTris(copy.nTris)
    {
        for (int i=0; i<copy.surfaces.size(); i++)
        {
            surfaces[i] = new surface(*copy.surfaces[i]);
        }
    }
    
    void readGeom(std::string geom_file);
    
    std::vector<surface*> getSurfaces() const {return surfaces;}
    
};


#endif /* defined(__CPanel__geometry__) */
