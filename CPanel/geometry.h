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
#include <fstream>
#include <array>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include "panelOctree.h"
#include "surface.h"
#include "liftingSurf.h"
#include "wake.h"
#include "bodyPanel.h"
#include "wakePanel.h"
#include "edge.h"
#include "inputParams.h"
#include "cpNode.h"

class geometry
{
    std::vector<liftingSurf*> liftingSurfs;
    std::vector<surface*> nonLiftingSurfs;
    std::vector<bodyPanel*> bPanels;
    std::vector<wakePanel*> wPanels;
    
    panelOctree pOctree;
    std::vector<cpNode*> nodes;
    std::vector<edge*> edges;
//    std::vector<cpNode*> TEnodes;
    short nNodes;
    short nTris;
    
    
    Eigen::MatrixXd B; // Source Influence Coefficient Matrix
    Eigen::MatrixXd A; // Doublet Influence Coefficient Matrix
    std::string infCoeffFile;

    void readTri(std::string tri_file, bool normFlag);
    std::vector<edge*> panEdges(const std::vector<cpNode*> &pNodes);
    edge* findEdge(cpNode* n1,cpNode* n2);
    void createSurfaces(const Eigen::MatrixXi &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID, std::vector<int> wakeIDs);
    void createOctree();
//    void setTEPanels();
//    void setTEnodes();
    void getLiftingSurfs(std::vector<surface*>& wakes, std::vector<surface*>& liftingSurfs);
    void setNeighbors(panel* p,int targetID);
//    void scanNode(panel* p, node<panel>* current, node<panel>* exception);
    bool isLiftingSurf(int currentID, std::vector<int> wakeIDs);
//    void correctWakeNodes(int wakeNodeStart);
    void correctWakeConnectivity(int wakeNodeStart,int wakeTriStart,Eigen::MatrixXi &connectivity);
    liftingSurf* getParentSurf(int wakeID);
    
    void setInfCoeff();
    Eigen::Vector4i interpIndices(std::vector<bodyPanel*> interpPans);
    
    bool infCoeffFileExists();
    void readInfCoeff();
    void writeInfCoeff();
    
public:
    geometry(inputParams* p)
    {
        std::stringstream temp;
        temp << p->geomFile->name << ".infCoeff";
        infCoeffFile = temp.str();
        readTri(p->geomFile->file, p->normFlag);
    }
    
    virtual ~geometry();
    
    geometry(const geometry& copy);
    
    geometry& operator=(const geometry &rhs);
    
    double pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d Vinf);
    
    short getNumberOfNodes() {return nNodes;}
    short getNumberOfTris() {return nTris;}
    std::vector<cpNode*> getNodes() {return nodes;}
    Eigen::MatrixXd getNodePnts();
    std::vector<liftingSurf*> getLiftingSurfs() {return liftingSurfs;}
    std::vector<surface*> getNonLiftingSurfs() {return nonLiftingSurfs;}
    std::vector<surface*> getSurfaces();
    panelOctree* getOctree() {return &pOctree;}
    std::vector<panel*> getPanels();
    std::vector<bodyPanel*>* getBodyPanels() {return &bPanels;}
    std::vector<wakePanel*>* getWakePanels() {return &wPanels;}
    std::vector<wake*> getWakes();
    Eigen::MatrixXd* getA() {return &A;}
    Eigen::MatrixXd* getB() {return &B;}
    
};


#endif /* defined(__CPanel__geometry__) */
