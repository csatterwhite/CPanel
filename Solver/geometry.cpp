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

geometry::geometry(std::string geom_file)
{
    readTri(geom_file);
}

void geometry::readTri(std::string tri_file)
{
    std::ifstream fid;
    fid.open(tri_file);
    if (fid.is_open())
    {
        fid >> nNodes >> nTris;
        Eigen::MatrixXi connectivity(nTris,3);
        Eigen::VectorXi allID(nTris);
        std::vector<int> surfIDs;
        std::vector<int> wakeIDs;
        std::vector<int> surfTypes;
        nodes.resize(nNodes,3);
        TEnodes.resize(nNodes,1);
        
        // Read XYZ Locations of Nodes
        for (int i=0; i<nNodes; i++)
        {
            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
        }
        
        findTEnodes();
        
        // Temporarily Store Connectivity
        for (int i=0; i<nTris; i++)
        {
            fid >> connectivity(i,0) >> connectivity(i,1) >> connectivity(i,2);
        }
        
        connectivity = connectivity.array()-1; //Adjust for 0 based indexing
        
        // Scan Surface IDs and collect Unique IDs
        for (int i=0; i<nTris; i++)
        {
            fid >> allID(i);
            if (i == 0 || allID(i) != allID(i-1))
            {
                if (allID(i) > 10000)
                {
                    wakeIDs.push_back(allID(i));
                }
                else
                {
                    surfIDs.push_back(allID(i));
                }
            }
        }
        createSurfaces(connectivity,allID,surfIDs,wakeIDs);
        createOctree();
        
        // Set neighbors
        
        std::vector<panel*> panels = getPanels();
        for (int i=0; i<panels.size(); i++)
        {
            panels[i]->setNeighbors(&pOctree);
        }
        
        // Set parents of trailing edge wake panels
        for (int i=0; i<liftingSurfs.size(); i++)
        {
            std::vector<wakePanel*> pans = liftingSurfs[i]->getWakePanels();
            for (int j=0; j<pans.size(); j++)
            {
                if (pans[j]->isTEpanel())
                {
                    pans[j]->setParentPanels();
                }
            }
        }
    }
    else
    {
        std::cout << "ERROR : File not found" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void geometry::findTEnodes()
{
    // Duplicate nodes are on TE because they are generated once for surface and once for wake. In some cases, there is floating point error that makes the points not exactly the same.  In these scenarios, the points are forced to be the exact same.
    Eigen::Vector3d vec;
    double eps = pow(10,-15);
    for (int i=0; i<nNodes; i++)
    {
        for (int j=0; j<nNodes; j++)
        {
            vec = nodes.row(i)-nodes.row(j);
            if (vec.norm()<eps && i != j)
            {
                nodes.row(j) = nodes.row(i);
                TEnodes(i) = true;
                TEnodes(j) = true;
            }
        }
    }
}

bool geometry::isLiftingSurf(int currentID, std::vector<int> wakeIDs)
{
    for (int i=0; i<wakeIDs.size(); i++)
    {
        if (wakeIDs[i]-10000 == currentID)
        {
            return true;
        }
    }
    return false;
}


void geometry::createSurfaces(Eigen::MatrixXi connectivity, Eigen::VectorXi allID, std::vector<int> surfIDs, std::vector<int> wakeIDs)
{
    surface* surf = nullptr;
    liftingSurf* surfL = nullptr;
    bool LS = false;
    for (int i=0; i<nTris; i++)
    {
        if (i==0 || allID(i)!=allID(i-1))
        {
            LS = isLiftingSurf(allID(i),wakeIDs);
            if (LS)
            {
                surfL = new liftingSurf(allID(i),&nodes);
                liftingSurfs.push_back(surfL);
            }
            else if (allID(i) <= 10000)
            {
                surf = new surface(allID(i),&nodes);
                nonLiftingSurfs.push_back(surf);
            }
            else
            {
                surfL = getParentSurf(allID(i));
            }
        }
        if (LS)
        {
            liftingSurfs.back()->addPanel(connectivity.row(i),TEnodes,allID(i));
        }
        else if (allID(i) <= 10000)
        {
            nonLiftingSurfs.back()->addPanel(connectivity.row(i),TEnodes);
        }
        else
        {
            surfL->addPanel(connectivity.row(i),TEnodes,allID(i));
        }
    }
}
         
liftingSurf* geometry::getParentSurf(int wakeID)
{
    for (int i=0; i<liftingSurfs.size(); i++)
    {
        if (liftingSurfs[i]->getID() == wakeID-10000)
        {
            return liftingSurfs[i];
        }
    }
    return nullptr;
}

void geometry::createOctree()
{
    std::vector<panel*> panels;
    std::vector<panel*> temp;
    std::vector<bodyPanel*> tempB;
    for (int i=0; i<liftingSurfs.size(); i++)
    {
        temp = liftingSurfs[i]->getAllPanels();
        panels.insert(panels.end(),temp.begin(),temp.end());
    }
    for (int i=0; i<nonLiftingSurfs.size(); i++)
    {
        tempB = nonLiftingSurfs[i]->getPanels();
        panels.insert(panels.end(),tempB.begin(),tempB.end());
    }
    pOctree.addData(panels);
}


std::vector<surface*> geometry::getSurfaces()
{
    std::vector<surface*> surfs;
    for (int i=0; i<nonLiftingSurfs.size(); i++)
    {
        surfs.push_back(nonLiftingSurfs[i]);
    }
    for (int i=0; i<liftingSurfs.size(); i++)
    {
        surfs.push_back(liftingSurfs[i]);
    }
    return surfs;
}

std::vector<panel*> geometry::getPanels()
{
    std::vector<panel*> panels;
    for (int i=0; i<liftingSurfs.size(); i++)
    {
        std::vector<panel*> temp = liftingSurfs[i]->getAllPanels();
        panels.insert(panels.end(),temp.begin(),temp.end());
    }
    for (int i=0; i<nonLiftingSurfs.size(); i++)
    {
        std::vector<bodyPanel*> temp = nonLiftingSurfs[i]->getPanels();
        panels.insert(panels.end(),temp.begin(),temp.end());
    }
    return panels;
}