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
#include "bodyPanel.h"
#include "wakePanel.h"
#include "edge.h"

class geometry
{
    double Sref;
    double bref;
    double cref;
    Eigen::Vector3d cg;
    
    std::vector<liftingSurf*> liftingSurfs;
    std::vector<surface*> nonLiftingSurfs;
    std::vector<bodyPanel*> bPanels;
    std::vector<wakePanel*> wPanels;
    
    panelOctree pOctree;
    Eigen::MatrixXd nodes;
    std::vector<edge*> edges;
    Eigen::Matrix<bool,Eigen::Dynamic,1> TEnodes;
    short nNodes;
    short nTris;
    
    
    Eigen::MatrixXd B; // Source Influence Coefficient Matrix
    Eigen::MatrixXd A; // Doublet Influence Coefficient Matrix

    void readTri(std::string tri_file, bool normFlag);
    std::vector<edge*> triEdges(const Eigen::VectorXi &indices);
    edge* findEdge(int i1,int i2);
    void createSurfaces(const Eigen::MatrixXi &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID, std::vector<int> wakeIDs);
    void createOctree();
    void setTEPanels();
    void getLiftingSurfs(std::vector<surface*>& wakes, std::vector<surface*>& liftingSurfs);
    void setNeighbors(panel* p,int targetID);
    void scanNode(panel* p, node<panel>* current, node<panel>* exception);
    bool isLiftingSurf(int currentID, std::vector<int> wakeIDs);
    void correctWakeNodes(int wakeNodeStart);
    void correctWakeConnectivity(int wakeNodeStart,int wakeTriStart,Eigen::MatrixXi &connectivity);
    liftingSurf* getParentSurf(int wakeID);
    
    void setInfCoeff();
    Eigen::Vector4i interpIndices(std::vector<bodyPanel*> interpPans);
    
public:
    geometry(std::string geomFile, bool normFlag, double Sref, double bref, double cref, Eigen::Vector3d cg) : Sref(Sref), bref(bref), cref(cref),cg(cg)
    {
        readTri(geomFile,normFlag);
    }
    
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
        for (int i=0; i<edges.size(); i++)
        {
            delete edges[i];
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
    
    double getSref() {return Sref;}
    double getbref() {return bref;}
    double getcref() {return cref;}
    
    Eigen::Vector3d getCG() {return cg;}
    short getNumberOfNodes() {return nNodes;}
    short getNumberOfTris() {return nTris;}
    Eigen::MatrixXd getNodes() {return nodes;}
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
