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
    
//    for (int i=0; i<nonLiftingSurfs.size(); i++)
//    {
//        delete nonLiftingSurfs[i];
//    }
//    nonLiftingSurfs.clear();
//    for (int i=0; i<liftingSurfs.size(); i++)
//    {
//        delete liftingSurfs[i];
//    }
//    liftingSurfs.clear();
    for (int i=0; i<surfaces.size(); i++)
    {
        delete surfaces[i];
    }
    surfaces.clear();
    for (int i=0; i<wakes.size(); i++)
    {
        delete wakes[i];
    }
    wakes.clear();
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
//    for (int i=0; i<copy.nonLiftingSurfs.size(); i++)
//    {
//        nonLiftingSurfs[i] = new surface(*copy.nonLiftingSurfs[i]);
//    }
//    for (int i=0; i<copy.liftingSurfs.size(); i++)
//    {
//        liftingSurfs[i] = new liftingSurf(*copy.liftingSurfs[i]);
//    }
    for (int i=0; i<copy.surfaces.size(); i++)
    {
        surfaces[i] = new surface(*copy.surfaces[i]);
    }
    for (int i=0; i<copy.wakes.size(); i++)
    {
        wakes[i] = new wake(*copy.wakes[i]);
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
//    for (int i=0; i<rhs.nonLiftingSurfs.size(); i++)
//    {
//        nonLiftingSurfs[i] = new surface(*rhs.nonLiftingSurfs[i]);
//    }
//    for (int i=0; i<rhs.liftingSurfs.size(); i++)
//    {
//        liftingSurfs[i] = new liftingSurf(*rhs.liftingSurfs[i]);
//    }
    for (int i=0; i<rhs.surfaces.size(); i++)
    {
        surfaces[i] = new surface(*rhs.surfaces[i]);
    }
    for (int i=0; i<rhs.wakes.size(); i++)
    {
        wakes[i] = new wake(*rhs.wakes[i]);
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
                if (allID(i) > 1000)
                {
                    wakeIDs.push_back(allID(i));
                }
                else
                {
                    surfIDs.push_back(allID(i));
                }
            }
            if (allID(i) > 1000 && allID(i-1) < 1000)
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
        
        std::cout << "\tNodes : " << nodes.size() << std::endl;
        std::cout << "\tEdges : " << edges.size() << std::endl;
        std::cout << "\tPanels : " << nTris << std::endl;
        
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
        
        
        if (wakes.size() > 1)
        {
            std::vector<wake*> newWakes;
            for (int i=0; i<wakes.size(); i++)
            {
                for (int j=i; j<wakes.size(); j++)
                {
                    if (wakes[i]->isSameWake(wakes[j]))
                    {
                        wakes[i]->mergeWake(wakes[j]);
                        delete wakes[j];
                        newWakes.push_back(wakes[i]);
                    }
                }
            }
            wakes = newWakes;
        }
        
        
        // Collect all panels in geometry
        
        std::vector<bodyPanel*> tempB;
        std::vector<wakePanel*> tempW;
        for (int i=0; i<surfaces.size(); i++)
        {
            tempB = surfaces[i]->getPanels();
            bPanels.insert(bPanels.begin(),tempB.begin(),tempB.end());
        }
        for (int i=0; i<wakes.size(); i++)
        {
            tempW = wakes[i]->getPanels();
            wPanels.insert(wPanels.begin(),tempW.begin(),tempW.end());
        }
        
        // Check panels for tip patches.  Needed to do 2D CHTLS to avoid nonphysical results near discontinuity at trailing edge.
        for (int i=0; i<bPanels.size(); i++)
        {
            bPanels[i]->setTipFlag();
        }
        for (int i=0; i<bPanels.size(); i++)
        {
            bPanels[i]->setCluster();
        }
        
//        std::cout << std::endl;
//        for (int i=0; i<bPanels.size(); i++)
//        {
//            if (bPanels[i]->isUpper())
//            {
//                std::cout << bPanels[i]->getCenter()(0) << "," << bPanels[i]->getCenter()(1) << "," << bPanels[i]->getCenter()(2) << ";" << std::endl;
//            }
//        }
//        std::cout << std::endl;
//        for (int i=0; i<bPanels.size(); i++)
//        {
//            if (bPanels[i]->isLower())
//            {
//                std::cout << bPanels[i]->getCenter()(0) << "," << bPanels[i]->getCenter()(1) << "," << bPanels[i]->getCenter()(2) << ";" << std::endl;
//            }
//        }
        
//        Eigen::Vector3d center;
//        std::vector<edge*> bodyPanEdges;
//        std::vector<edge*> wakePanEdges;
//        for (int i=0; i<bPanels.size(); i++)
//        {
//            center = bPanels[i]->getCenter();
//            if (center(0) > 5.45 && center(0) < 5.74 && center(1) > -0.77 && center(1) < -0.65 && center(2) > -0.59 && center(2) < -0.42)
//            {
//                bodyPanEdges = bPanels[i]->getEdges();
//            }
//        }
//        
//        for (int i=0; i<wPanels.size(); i++)
//        {
//            center = wPanels[i]->getCenter();
//            if (center(0) > 5.45 && center(0) < 5.74 && center(1) > -0.92 && center(1) < -0.65 && center(2) > -0.59 && center(2) < -0.58)
//            {
//                wakePanEdges = wPanels[i]->getEdges();
//            }
//        }
        
        // Calculate influence coefficient matrices
        
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
        std::cout << "ERROR : Geometry file not found" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void geometry::correctWakeConnectivity(int wakeNodeStart,int wakeTriStart,Eigen::MatrixXi &connectivity)
{
    Eigen::Vector3d vec;
    Eigen::Matrix<double,Eigen::Dynamic,2> indReps; // [toReplace, replaceWith]
    double tol = shortestEdge(connectivity);
    int count = 0;
    for (int i=0; i<wakeNodeStart; i++)
    {
        for (int j=wakeNodeStart; j<nNodes; j++)
        {
            vec = nodes[i]->getPnt()-nodes[j]->getPnt();
            if (vec.lpNorm<Eigen::Infinity>() < tol)
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

double geometry::shortestEdge(const Eigen::MatrixXi &connectivity)
{
    int nrows,ncols;
    Eigen::Vector3d p1,p2;
    nrows = (int)connectivity.rows();
    ncols = (int)connectivity.cols();
    double l;
    double shortest = 100000;
    for (int i=0; i<nrows; i++)
    {
        for (int j=0; j<ncols; j++)
        {
            p1 = nodes[connectivity(i,j)]->getPnt();
            if (j < ncols-1)
            {
                p2 = nodes[connectivity(i,j+1)]->getPnt();
            }
            else
            {
                p2 = nodes[connectivity(i,0)]->getPnt();
            }
            l = (p2-p1).norm();
            if (l < shortest)
            {
                shortest = l;
            }
        }
    }
    return shortest;
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
    surface* s = nullptr;
    wake* w = nullptr;
    bodyPanel* bPan;
    wakePanel* wPan;
    std::vector<edge*> pEdges;
    for (int i=0; i<nTris; i++)
    {
        std::vector<cpNode*> pNodes;
        pNodes.push_back(nodes[connectivity(i,0)]);
        pNodes.push_back(nodes[connectivity(i,1)]);
        pNodes.push_back(nodes[connectivity(i,2)]);
        pEdges = panEdges(pNodes); //Create edge or find edge that already exists
        if (allID(i) <= 1000)
        {
            if (i==0 || allID(i) != allID(i-1))
            {
                s = new surface(allID(i),this);
                surfaces.push_back(s);
            }
            bPan = new bodyPanel(pNodes,pEdges,norms.row(i),s,allID(i));
            s->addPanel(bPan);
        }
        else
        {
            if (i==0 || allID(i)!=allID(i-1))
            {
                w = new wake(allID(i));
                wakes.push_back(w);
            }
            wPan = new wakePanel(pNodes,pEdges,norms.row(i),w,allID(i));
            w->addPanel(wPan);
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

void geometry::createOctree()
{
    std::vector<panel*> panels;
    std::vector<wakePanel*> tempW;
    std::vector<bodyPanel*> tempB;

    for (int i=0; i<surfaces.size(); i++)
    {
        tempB = surfaces[i]->getPanels();
        panels.insert(panels.end(),tempB.begin(),tempB.end());
    }
    for (int i=0; i<wakes.size(); i++)
    {
        tempW = wakes[i]->getPanels();
        panels.insert(panels.end(),tempW.begin(),tempW.end());
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

    if (writeCoeffFlag)
    {
        writeInfCoeff();
    }
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
    return surfaces;
}

std::vector<wake*> geometry::getWakes()
{
    return wakes;
}

std::vector<panel*> geometry::getPanels()
{
    std::vector<panel*> panels;

    for (int i=0; i<surfaces.size(); i++)
    {
        std::vector<bodyPanel*> temp = surfaces[i]->getPanels();
        panels.insert(panels.end(),temp.begin(),temp.end());
    }
    for (int i=0; i<wakes.size(); i++)
    {
        std::vector<wakePanel*> temp = wakes[i]->getPanels();
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

Eigen::Vector3d geometry::pntVelocity(const Eigen::Vector3d &pnt, const Eigen::Vector3d &Vinf, double PG)
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
    
    vel(0) /= PG; // Prandtl-Glauert Correction
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
