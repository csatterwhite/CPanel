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
        
        // Temporarily Store Connectivity
        for (int i=0; i<nTris; i++)
        {
            fid >> connectivity(i,0) >> connectivity(i,1) >> connectivity(i,2);
        }
        
        connectivity = connectivity.array()-1; //Adjust for 0 based indexing
        
        // Scan Surface IDs and collect Unique IDs
        int wakeNodeStart = nNodes;
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
            if (allID(i) > 10000 && allID(i-1) < 10000)
            {
                wakeNodeStart = connectivity.row(i).minCoeff();
            }
        }
        createSurfaces(connectivity,allID,surfIDs,wakeIDs);
        createOctree();
        if (wakeIDs.size() > 0)
        {
            correctWakeNodes(wakeNodeStart);
        }
        
        // Set neighbors
        for (int i=0; i<liftingSurfs.size(); i++)
        {
            // Set wake neighbors first to deal with special cases
            liftingSurfs[i]->getWake()->setNeighbors(&pOctree);
        }
        
        std::vector<surface*> surfs = getSurfaces();
        for (int i=0; i<surfs.size(); i++)
        {
            surfs[i]->setNeighbors(&pOctree);
        }
    }
    else
    {
        std::cout << "ERROR : File not found" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void geometry::correctWakeNodes(int wakeNodeStart)
{
    // Duplicate nodes are on TE because they are generated once for surface and once for wake. In some cases, there is error that makes the points not exactly the same.  In these scenarios, the points are forced to be the exact same. The highest error scene has been on the order of 10^-5.  This error comes from VSPs node generation and is not part of CPanel
    Eigen::Vector3d vec;
    double diff = pow(10,-3);
    for (int i=0; i<wakeNodeStart; i++)
    {
        for (int j=wakeNodeStart; j<nNodes; j++)
        {
            vec = nodes.row(i)-nodes.row(j);
            if (vec.lpNorm<Eigen::Infinity>() < diff)
            {
                nodes.row(j) = nodes.row(i);
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

std::vector<wakePanel*> geometry::getWakePanels()
{
    std::vector<wakePanel*> wps;
    for (int i=0; i<liftingSurfs.size(); i++)
    {
        wps.insert(wps.end(),liftingSurfs[i]->getWakePanels().begin(),liftingSurfs[i]->getWakePanels().end());
    }
    return wps;
}