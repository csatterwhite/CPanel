//
//  geometry.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "geometry.h"

// Destructor //

geometry::~geometry()
{
    
    for (int i=0; i<nonLiftingSurfs.size(); i++)
    {
        delete nonLiftingSurfs[i];
    }
    nonLiftingSurfs.clear();
    for (int i=0; i<liftingSurfs.size(); i++)
    {
        delete liftingSurfs[i];
    }
    liftingSurfs.clear();
    for (int i=0; i<edges.size(); i++)
    {
        delete edges[i];
    }
    edges.clear();
    for (int i=0; i<nodes.size(); i++)
    {
        delete nodes[i];
    }
    nodes.clear();
//    for (int i=0; i<TEnodes.size(); i++)
//    {
//        delete TEnodes[i];
//    }
//    TEnodes.clear();
    
}

// Copy Constructor //

geometry::geometry(const geometry& copy) : pOctree(copy.pOctree), nNodes(copy.nNodes), nTris(copy.nTris), A(copy.A), B(copy.B), infCoeffFile(copy.infCoeffFile)
{
    for (int i=0; i<copy.nonLiftingSurfs.size(); i++)
    {
        nonLiftingSurfs[i] = new surface(*copy.nonLiftingSurfs[i]);
    }
    for (int i=0; i<copy.liftingSurfs.size(); i++)
    {
        liftingSurfs[i] = new liftingSurf(*copy.liftingSurfs[i]);
    }
    for (int i=0; i<copy.bPanels.size(); i++)
    {
        bPanels[i] = new bodyPanel(*copy.bPanels[i]);
    }
    for (int i=0; i<copy.wPanels.size(); i++)
    {
        wPanels[i] = new wakePanel(*copy.wPanels[i]);
    }
    for (int i=0; i<copy.nodes.size(); i++)
    {
        nodes[i] = new cpNode(*copy.nodes[i]);
    }
    for (int i=0; i<copy.edges.size(); i++)
    {
        edges[i] = new edge(*copy.edges[i]);
    }
//    for (int i=0; i<copy.TEnodes.size(); i++)
//    {
//        TEnodes[i] = new cpNode(*copy.TEnodes[i]);
//    }
}

// Assignment Operator //

geometry& geometry::operator=(const geometry &rhs)
{
    if (this == &rhs)
    {
        return (*this);
    }
    
    pOctree = rhs.pOctree;
    nNodes = rhs.nNodes;
    nTris = rhs.nTris;
    A = rhs.A;
    B = rhs.B;
    infCoeffFile = rhs.infCoeffFile;
    
    // Deep Copy of pointers
    for (int i=0; i<rhs.nonLiftingSurfs.size(); i++)
    {
        nonLiftingSurfs[i] = new surface(*rhs.nonLiftingSurfs[i]);
    }
    for (int i=0; i<rhs.liftingSurfs.size(); i++)
    {
        liftingSurfs[i] = new liftingSurf(*rhs.liftingSurfs[i]);
    }
    for (int i=0; i<rhs.bPanels.size(); i++)
    {
        bPanels[i] = new bodyPanel(*rhs.bPanels[i]);
    }
    for (int i=0; i<rhs.wPanels.size(); i++)
    {
        wPanels[i] = new wakePanel(*rhs.wPanels[i]);
    }
    for (int i=0; i<rhs.nodes.size(); i++)
    {
        nodes[i] = new cpNode(*rhs.nodes[i]);
    }
    for (int i=0; i<rhs.edges.size(); i++)
    {
        edges[i] = new edge(*rhs.edges[i]);
    }
    
    return *this;
}

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
        
        // Read XYZ Locations of Nodes
        Eigen::Vector3d pnt;
        cpNode* n;
        for (int i=0; i<nNodes; i++)
        {
            fid >> pnt(0) >> pnt(1) >> pnt(2);
            n = new cpNode(pnt,i);
            nodes.push_back(n);
        }
        
        // Temporarily Store Connectivity
        for (int i=0; i<nTris; i++)
        {
            fid >> connectivity(i,0) >> connectivity(i,1) >> connectivity(i,2);
        }
        
        connectivity = connectivity.array()-1; //Adjust for 0 based indexing
        
        // Scan Surface IDs and collect Unique IDs
        int wakeNodeStart = nNodes;
        int wakeTriStart = nTris;
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
                wakeTriStart = i;
            }
        }
        
        if (wakeIDs.size() > 0)
        {
            correctWakeConnectivity(wakeNodeStart, wakeTriStart, connectivity);
        }
        
        // Read in Normals if included in input file
        Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(nTris,3);
        if (normFlag)
        {
            std::cout << "Reading Bezier Normals from Geometry File..." << std::endl;
            for (int i=0; i<nTris; i++)
            {
                fid >> norms(i,0) >> norms(i,1) >> norms(i,2);
            }
        }
        
        std::cout << "Generating Panel Geometry..." << std::endl;
        
        createSurfaces(connectivity,norms,allID,wakeIDs);
        
        for (int i=0; i<edges.size(); i++)
        {
            if (edges[i]->getBodyPans().size() == 0 && edges[i]->isTE())
            {
                assert(edges[i]->isTE() != true);
            }
        }
        
        // Erase duplicate node pointers
        std::sort( nodes.begin(), nodes.end() );
        nodes.erase( std::unique( nodes.begin(), nodes.end() ), nodes.end() );
        
        for (int i=0; i<nodes.size(); i++)
        {
            nodes[i]->setIndex(i);
        }
        
        std::cout << "Building Octree..." << std::endl;

        createOctree();
        
        // Set neighbors
        std::cout << "Finding Panel Neighbors..." << std::endl;
        
        for (int i=0; i<edges.size(); i++)
        {
            edges[i]->setNeighbors();
        }
        for (int i=0; i<wPanels.size(); i++)
        {
            if (wPanels[i]->isTEpanel())
            {
                wPanels[i]->setParentPanels();
            }
        }
        
//        setTEnodes();
        
        bool read = false;
        
        if (infCoeffFileExists())
        {
            std::string in;
            std::cout << "\nInfluence Coefficients have already been calculated for a geometry with this name, would you like to use these coefficients?" << std::endl;
            std::cout << "\t< Y > - Yes, use coefficients." << std::endl;
            std::cout << "\t< N > - No, recalculate them." << std::endl;
            std::cin >> in;
            std::cout << std::endl;
            if (in == "Y")
            {
                readInfCoeff();
                read = true;
            }
        }
        
        if (!read)
        {
            std::cout << "Building Influence Coefficient Matrix..." << std::endl;
            setInfCoeff();
        }
    }
    else
    {
        std::cout << "ERROR : File not found" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void geometry::correctWakeConnectivity(int wakeNodeStart,int wakeTriStart,Eigen::MatrixXi &connectivity)
{
    Eigen::Vector3d vec;
    Eigen::Matrix<double,Eigen::Dynamic,2> indReps; // [toReplace, replaceWith]
    double diff = pow(10,-3);
    int count = 0;
    for (int i=0; i<wakeNodeStart; i++)
    {
        for (int j=wakeNodeStart; j<nNodes; j++)
        {
            vec = nodes[i]->getPnt()-nodes[j]->getPnt();
            if (vec.lpNorm<Eigen::Infinity>() < diff)
            {
                nodes[j] = nodes[i];
                count++;
                indReps.conservativeResize(count,2);
                indReps(count-1,0) = j;
                indReps(count-1,1) = i;
            }
        }
    }
    
    for (int i=wakeTriStart; i<nTris; i++)
    {
        for (int j=0; j<connectivity.cols(); j++)
        {
            for (int k=0; k<indReps.rows(); k++)
            {
                if (connectivity(i,j) == indReps(k,0))
                {
                    connectivity(i,j) = indReps(k,1);
                }
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
    std::vector<edge*> pEdges;
    bool LS = false;
    for (int i=0; i<nTris; i++)
    {
        std::vector<cpNode*> pNodes;
        pNodes.push_back(nodes[connectivity(i,0)]);
        pNodes.push_back(nodes[connectivity(i,1)]);
        pNodes.push_back(nodes[connectivity(i,2)]);
        pEdges = panEdges(pNodes);
        if (i==0 || allID(i)!=allID(i-1))
        {
            LS = isLiftingSurf(allID(i),wakeIDs);
            if (LS)
            {
                surfL = new liftingSurf(allID(i),this);
                liftingSurfs.push_back(surfL);
            }
            else if (allID(i) <= 10000)
            {
                surf = new surface(allID(i),this);
                nonLiftingSurfs.push_back(surf);
            }
            else
            {
                surfL = getParentSurf(allID(i));
            }
        }
        if (LS)
        {
            bPan = new bodyPanel(pNodes,pEdges,norms.row(i),liftingSurfs.back(),allID(i),true);
            liftingSurfs.back()->addPanel(bPan);
            bPanels.push_back(bPan);
        }
        else if (allID(i) <= 10000)
        {
            bPan = new bodyPanel(pNodes,pEdges,norms.row(i),nonLiftingSurfs.back(),allID(i),false);
            nonLiftingSurfs.back()->addPanel(bPan);
            bPanels.push_back(bPan);
        }
        else
        {
            wPan = new wakePanel(pNodes,pEdges,norms.row(i),allID(i),surfL->getWake());
            surfL->addPanel(wPan);
            wPanels.push_back(wPan);
        }
    }
}

std::vector<edge*> geometry::panEdges(const std::vector<cpNode*>  &pNodes)
{
    int i1,i2;
    std::vector<edge*> triEdges;
    edge* e;
    for (int i=0; i<pNodes.size(); i++)
    {
        i1 = i;
        if (i == pNodes.size()-1)
        {
            i2 = 0;
        }
        else
        {
            i2 = i+1;
        }
        e = findEdge(pNodes[i1],pNodes[i2]);
        triEdges.push_back(e);
    }
    return triEdges;
}

edge* geometry::findEdge(cpNode* n1,cpNode* n2)
{
    for (int i=0; i<edges.size(); i++)
    {
        if (edges[i]->sameEdge(n1, n2))
        {
            return edges[i];
        }
    }
    
    // If edge doesn't exist, create one
    edge* e = new edge(n1,n2,this);
    edges.push_back(e);
    return e;
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
    
    A.resize(nBodyPans,nBodyPans);
    B.resize(nBodyPans,nBodyPans);
    
    Eigen::VectorXi percentage(9);
    percentage << 10,20,30,40,50,60,70,80,90;
    
    for (int j=0; j<nBodyPans; j++)
    {
        for (int i=0; i<nBodyPans; i++)
        {
            bPanels[j]->panelPhiInf(bPanels[i]->getCenter(),B(i,j),A(i,j));
        }
        for (int i=0; i<percentage.size(); i++)
        {
            if ((100*j/nPans) <= percentage(i) && 100*(j+1)/nPans > percentage(i))
            {
                std::cout << percentage(i) << "%\t" << std::flush;
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
                std::cout << percentage(i) << "%\t" << std::flush;
            }
        }
        
    }
    std::cout << "Complete" << std::endl;
    writeInfCoeff();
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

bool geometry::infCoeffFileExists()
{
    boost::filesystem::path p = boost::filesystem::current_path().string()+"/" + infCoeffFile;
    if (boost::filesystem::exists(p))
    {
        std::ifstream fid;
        fid.open(p.string());
        int nPans;
        fid >> nPans;
        fid.close();
        if (nPans != bPanels.size())
        {
            return false;
        }
        return true;
    }
    
    return false;
}

void geometry::readInfCoeff()
{
    std::cout << "Reading Influence Coefficients from " << infCoeffFile << "..." << std::endl;
    
    std::ifstream fid;
    fid.open(infCoeffFile);
    int nPans;
    fid >> nPans;
    A.resize(nPans,nPans);
    B.resize(nPans,nPans);
    for (int i=0; i<bPanels.size(); i++)
    {
        for (int j=0; j<bPanels.size(); j++)
        {
            fid >> A(i,j);
        }
    }
    for (int i=0; i<bPanels.size(); i++)
    {
        for (int j=0; j<bPanels.size(); j++)
        {
            fid >> B(i,j);
        }
    }
    fid.close();
}

void geometry::writeInfCoeff()
{
    std::cout << "Writing Influence Coefficients to " << infCoeffFile << "..." << std::endl;
    std::ofstream fid;
    fid.open(infCoeffFile);
    fid << bPanels.size() << "\n";
    for (int i=0; i<bPanels.size(); i++)
    {
        for (int j=0; j<bPanels.size(); j++)
        {
            fid << A(i,j) << "\t";
        }
        fid << "\n";
    }
    for (int i=0; i<bPanels.size(); i++)
    {
        for (int j=0; j<bPanels.size(); j++)
        {
            fid << B(i,j) << "\t";
        }
        fid << "\n";
    }
    fid.close();
}

Eigen::MatrixXd geometry::getNodePnts()
{
    Eigen::MatrixXd nodePnts(nodes.size(),3);
    for (int i=0; i<nodes.size(); i++)
    {
        nodePnts.row(i) = nodes[i]->getPnt();
    }
    return nodePnts;
}

double geometry::pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf)
{
    double pot = 0;
    for (int i=0; i<bPanels.size(); i++)
    {
        pot += bPanels[i]->panelPhi(pnt);
    }
    for (int i=0; i<wPanels.size(); i++)
    {
        pot += wPanels[i]->panelPhi(pnt);
    }
    pot += Vinf.dot(pnt);
    return pot;
}

Eigen::Vector3d geometry::pntVelocity(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf)
{
    Eigen::Vector3d vel = Eigen::Vector3d::Zero();
    for (int i=0; i<bPanels.size(); i++)
    {
        vel += bPanels[i]->panelV(pnt);
    }
    for (int i=0; i<wPanels.size(); i++)
    {
        vel += wPanels[i]->panelV(pnt);
    }
    vel += Vinf;
    return vel;
}


void geometry::clusterCheck()
{
    std::ofstream fid;
    fid.open("ClusterCheck.txt");
    int index;
    std::vector<bodyPanel*> clust;
    for (int i=0; i<bPanels.size(); i++)
    {
        clust = bPanels[i]->getCluster();
        for (int j=0; j<clust.size(); j++)
        {
            index = (int)std::distance(bPanels.begin(),std::find(bPanels.begin(), bPanels.end(), clust[j]));
            fid << index+1 << "\t";
        }
        fid << "\n";
    }
    fid.close();
}
