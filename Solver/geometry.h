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
#include "panelOctree.h"
#include "surface.h"
#include "liftingSurf.h"
#include "wake.h"
#include "wakePanel.h"
#include "bodyPanel.h"

class geometry
{
    std::vector<liftingSurf*> liftingSurfs;
    std::vector<surface*> nonLiftingSurfs;
    panelOctree pOctree;
    Eigen::MatrixXd nodes;
    Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes;
    short nNodes;
    short nTris;

    void readTri(std::string tri_file, bool normFlag);
    void createSurfaces(const Eigen::MatrixXi &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID, std::vector<int> wakeIDs);
    void createOctree();
    void setTEPanels();
    void getLiftingSurfs(std::vector<surface*>& wakes, std::vector<surface*>& liftingSurfs);
    void setNeighbors(panel* p,int targetID);
    void scanNode(panel* p, node<panel>* current, node<panel>* exception);
    bool isLiftingSurf(int currentID, std::vector<int> wakeIDs);
    void correctWakeNodes(int wakeNodeStart);
    liftingSurf* getParentSurf(int wakeID);
    
public:
    geometry(std::string tri_file, bool normFlag);
    
    ~geometry()
    {
        for (int i=0; i<nonLiftingSurfs.size(); i++)
        {
            delete nonLiftingSurfs[i];
        }
        for (int i=0; i<liftingSurfs.size(); i++)
        {
            delete liftingSurfs[i];
        }
    }
    
    geometry(const geometry& copy) : pOctree(copy.pOctree), nodes(copy.nodes), nNodes(copy.nNodes), nTris(copy.nTris)
    {
        for (int i=0; i<copy.nonLiftingSurfs.size(); i++)
        {
            nonLiftingSurfs[i] = new surface(*copy.nonLiftingSurfs[i]);
        }
        for (int i=0; i<copy.liftingSurfs.size(); i++)
        {
            liftingSurfs[i] = new liftingSurf(*copy.liftingSurfs[i]);
        }
    }
    
    short getNumberOfNodes() {return nNodes;}
    short getNumberOfTris() {return nTris;}
    Eigen::MatrixXd getNodes() {return nodes;}
    std::vector<liftingSurf*> getLiftingSurfs() {return liftingSurfs;}
    std::vector<surface*> getNonLiftingSurfs() {return nonLiftingSurfs;}
    std::vector<surface*> getSurfaces();
    panelOctree* getOctree() {return &pOctree;}
    std::vector<panel*> getPanels();
    std::vector<wakePanel*> getWakePanels();
    std::vector<wake*> getWakes();
    
};


#endif /* defined(__CPanel__geometry__) */
