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
        
        connectivity = connectivity.array()-1; //Adjust for 0 based indexing
        
        fixDuplicateNodes(connectivity);
        
        // Scan Surface IDs and collect Unique IDs
        for (int i=0; i<nTris; i++)
        {
            fid >> surfID(i);
        }
        createSurfaces(connectivity,surfID);
        createOctree();
        setTEPanels();
        
        short count = 0;
        for (int i=0; i<surfaces.size(); i++)
        {
            std::vector<panel*> temp = surfaces[i]->getPanels();
            for (int j=0; j<temp.size(); j++)
            {
                if (temp[j]->isTEpanel())
                {
                    count++;
                }
            }
        }
        std::cout << count << std::endl;
    }
}

void geometry::fixDuplicateNodes(Eigen::MatrixXi &connectivity)
{
    for (int i=0; i<nodes.rows()-1; i++)
    {
        for (int j=i+1; j<nodes.rows(); j++)
        {
            if (isSameNode(nodes.row(i),nodes.row(j)))
            {
                changeIndex(connectivity, j, i);
            }
        }
    }
}

bool geometry::isSameNode(Eigen::Vector3d p1, Eigen::Vector3d p2)
{
    if (p1(0) == p2(0) && p1(1) == p2(1) && p1(2) == p2(2))
    {
        return true;
    }
    else
    {
        return false;
    }
}

void geometry::changeIndex(Eigen::MatrixXi &connectivity, int toReplace, int replaceWith)
{
    for (int i=0; i<connectivity.rows(); i++)
    {
        for (int j=0; j<3; j++)
        {
            if (connectivity(i,j) == toReplace)
            {
                connectivity(i,j) = replaceWith;
            }
        }
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

void geometry::setTEPanels()
{
    std::vector<surface*> wakes;
    std::vector<surface*> liftingSurfs;
    getLiftingSurfs(wakes,liftingSurfs);
    for (int i=0; i<wakes.size(); i++)
    {
        std::vector<panel*> wakePanels = wakes[i]->getPanels();
        short targetID = liftingSurfs[i]->getID();
        for (int j=0; j<wakePanels.size(); j++)
        {
            setNeighbors(wakePanels[j],targetID);
        }
    }
}

void geometry::getLiftingSurfs(std::vector<surface*>& wakes, std::vector<surface*>& liftingSurfs)
{
    for (int i=0; i<surfaces.size(); i++)
    {
        if (surfaces[i]->getID()>10000)
        {
            wakes.push_back(surfaces[i]);
        }
    }
    for (int i=0; i<wakes.size(); i++)
    {
        for (int j=0; j<surfaces.size(); j++)
        {
            if (surfaces[j]->getID() == wakes[i]->getID()-10000)
            {
                liftingSurfs.push_back(surfaces[j]);
            }
        }
    }
}

void geometry::setNeighbors(panel* p, int targetID)
{
    node<panel>* currentNode = pOctree.findNodeContainingMember(p);
    scanNode(p,currentNode,NULL);
    bool flag = false; //Flags panels that are known to not have any trailing edge panels for neighbors;
    while (currentNode != pOctree.getRootNode())
    {
        scanNode(p,currentNode->getParent(),currentNode);
        currentNode = currentNode->getParent();
        if (p->getNeighbors().size() == 3)
        {
            std::vector<panel*> neighbors = p->getNeighbors();
            short count = 0; //Count neighbors that are not on lifting surface.
            for (int i=0; i<3; i++)
            {
                if (neighbors[i]->getID() != targetID)
                {
                    count++;
                }
            }
            if (count == 3)
            {
                flag = true;
                break; //If there are three neighbors and none of them are on the lifting surface, the panel is not touching the trailing edge.
            }
        }
    }

    if (!flag && p->getNeighbors().size()>=3)
    {
        std::vector<panel*> neighbors = p->getNeighbors();
        for (int i=0; i<neighbors.size(); i++)
        {
            if (neighbors[i]->getID() == targetID)
            {
                neighbors[i]->setTE();
            }
        }
    }
    
}

void geometry::scanNode(panel* p, node<panel>* current, node<panel>* exception)
{
    std::vector<panel*> nodeMembers = current->getMembers(exception);
    for (int i=0; i<nodeMembers.size(); i++)
    {
        p->checkNeighbor(nodeMembers[i]);
    }
}
