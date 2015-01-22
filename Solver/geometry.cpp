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

void geometry::readTri(std::string tri_file, bool normFlag)
{
    std::ifstream fid;
    fid.open(tri_file);
    if (fid.is_open())
    {
        std::cout << "Reading Geometry..." << std::endl;
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
        
        if (wakeIDs.size() > 0)
        {
            correctWakeNodes(wakeNodeStart);
        }
        
        // Read in Normals if included in input file
        Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(nTris,3);
        if (normFlag)
        {
            for (int i=0; i<nTris; i++)
            {
                fid >> norms(i,0) >> norms(i,1) >> norms(i,2);
            }
        }
        
        std::cout << "Building Octree..." << std::endl;
        
        createSurfaces(connectivity,norms,allID,wakeIDs);
        createOctree();
        
        // Set neighbors
        for (int i=0; i<liftingSurfs.size(); i++)
        {
            // Set wake neighbors first to deal with special cases
            liftingSurfs[i]->getWake()->setNeighbors(&pOctree);
        }
        
        std::cout << "Finding Panel Neighbors..." << std::endl;
        
        std::vector<surface*> surfs = getSurfaces();
        for (int i=0; i<surfs.size(); i++)
        {
            surfs[i]->setNeighbors(&pOctree);
        }
        
        std::cout << "Building Influence Coefficient Matrix..." << std::endl;
        setInfCoeff();
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


void geometry::createSurfaces(const Eigen::MatrixXi &connectivity, const Eigen::MatrixXd &norms, const Eigen::VectorXi &allID, std::vector<int> wakeIDs)
{
    surface* surf = nullptr;
    liftingSurf* surfL = nullptr;
    bodyPanel* bPan;
    wakePanel* wPan;
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
            bPan = new bodyPanel(connectivity.row(i),&nodes,norms.row(i),allID(i),true);
            liftingSurfs.back()->addPanel(bPan);
            bPanels.push_back(bPan);
        }
        else if (allID(i) <= 10000)
        {
            bPan = new bodyPanel(connectivity.row(i),&nodes,norms.row(i),allID(i),false);
            nonLiftingSurfs.back()->addPanel(bPan);
            bPanels.push_back(bPan);
        }
        else
        {
            wPan = new wakePanel(connectivity.row(i),&nodes,norms.row(i),allID(i),surfL->getWake());
            surfL->addPanel(wPan);
            wPanels.push_back(wPan);
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

void geometry::setInfCoeff()
{
    // Construct doublet and source influence coefficient matrices for body panels
    int nBodyPans = (int)bPanels.size();
    int nWakePans = (int)wPanels.size();
    int nPans = nBodyPans+nWakePans;
    Eigen::VectorXd sigmas(nBodyPans);
    
    A.resize(nBodyPans,nBodyPans);
    B.resize(nBodyPans,nBodyPans);
    
    Eigen::VectorXi percentage(9);
    percentage << 10,20,30,40,50,60,70,80,90;
    
    for (int j=0; j<nBodyPans; j++)
    {
        for (int i=0; i<nBodyPans; i++)
        {
            bPanels[j]->panelPhiInf(bPanels[i]->getCenter(),B(i,j),A(i,j));
            if (i==j)
            {
                A(i,j) = -0.5;
            }
            if (j==0)
            {
                sigmas(i) = bPanels[i]->getSigma();
            }
        }
        for (int i=0; i<percentage.size(); i++)
        {
            if ((100*j/nPans) <= percentage(i) && 100*(j+1)/nPans > percentage(i))
            {
                std::cout << percentage(i) << "%" << std::endl;
            }
        }
    }
    
    for (int i=0; i<nBodyPans; i++)
    {
        bPanels[i]->setIndex(i);
    }
    
    std::vector<bodyPanel*> interpPans(4); // [Upper1 Lower1 Upper2 Lower2]  Panels that start the bounding wakelines of the wake panel.  Doublet strength is constant along wakelines (muUpper-muLower) and so the doublet strength used for influence of wake panel is interpolated between wakelines.
    double interpCoeff;
    double influence;
    Eigen::Vector4i indices;
    
    for (int j=0; j<nWakePans; j++)
    {
        wPanels[j]->interpPanels(interpPans,interpCoeff);
        indices = interpIndices(interpPans);
        for (int i=0; i<nBodyPans; i++)
        {
            influence = wPanels[j]->dubPhiInf(bPanels[i]->getCenter());
            A(i,indices(0)) += influence*(1-interpCoeff);
            A(i,indices(1)) += influence*(interpCoeff-1);
            A(i,indices(2)) += influence*interpCoeff;
            A(i,indices(3)) -= influence*interpCoeff;
        }
        for (int i=0; i<percentage.size(); i++)
        {
            if ((100*(nBodyPans+j)/nPans) <= percentage(i) && 100*(nBodyPans+j+1)/nPans > percentage(i))
            {
                std::cout << percentage(i) << "%" << std::endl;
            }
        }
        
    }
    std::cout << "Complete" << std::endl;
}

Eigen::Vector4i geometry::interpIndices(std::vector<bodyPanel*> interpPans)
{
    Eigen::Vector4i indices;
    for (int i=0; i<interpPans.size(); i++)
    {
        indices(i) = interpPans[i]->getIndex();
    }
    return indices;
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

std::vector<wake*> geometry::getWakes()
{
    std::vector<wake*> wakes;
    for (int i=0; i<liftingSurfs.size(); i++)
    {
        wakes.push_back(liftingSurfs[i]->getWake());
    }
    return wakes;
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