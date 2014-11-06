//
//  runCase.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "runCase.h"

runCase::runCase(geometry *geom,double V,double alpha,double beta,std::string outFile) : geom(geom), outFile(outFile)
{
    Vinf = windToBody(V,alpha,beta);
    setSourceStrengths();
    solveMatrixEq();
}

Eigen::Vector3d runCase::windToBody(double V, double alpha, double beta)
{
    alpha *= M_PI/180;
    beta *= M_PI/180;
    
    Eigen::Matrix3d T;
    T(0,0) = cos(alpha)*cos(beta);
    T(0,1) = cos(alpha)*sin(beta);
    T(0,2) = -sin(alpha);
    T(1,0) = -sin(beta);
    T(1,1) = cos(beta);
    T(1,2) = 0;
    T(2,0) = sin(alpha)*cos(beta);
    T(2,1) = sin(alpha)*sin(beta);
    T(2,2) = cos(alpha);
    Eigen::Vector3d Vt;
    Eigen::Vector3d Vel;
    Vt << V,0,0;

    Vel = T*Vt;
    
    return Vel;
}

void runCase::setSourceStrengths()
{
    std::vector<surface*> surfaces = geom->getSurfaces();
    std::vector<bodyPanel*> panels;
    for (int i=0; i<surfaces.size(); i++)
    {
        panels = surfaces[i]->getPanels();
        for (int j=0; j<panels.size(); j++)
        {
            panels[j]->setSigma(Vinf,0);
//            std::cout << panels[j]->getSigma() << std::endl;
//            std::cout << panels[j]->getCenter() << std::endl;
        }
    }
}

void runCase::solveMatrixEq()
{
    std::vector<surface*> surfaces = geom->getSurfaces();
    std::vector<bodyPanel*> tempB;
    std::vector<wakePanel*> tempW;
    
    for (int i=0; i<surfaces.size(); i++)
    {
        tempB = surfaces[i]->getPanels();
        for (int j=0; j<tempB.size(); j++)
        {
            if (tempB[j]->isTEpanel())
            {
                // Insert TEpanels at the beginning of the vector to speed up search for index location of TEpanel shedding a wakeline needed for computing the influence of a wake panel.
                bPanels.insert(bPanels.begin(), tempB[j]);
            }
            else
            {
                bPanels.push_back(tempB[j]);
            }
        }
    }
    
    for (int i=0; i<geom->getLiftingSurfs().size(); i++)
    {
        tempW = geom->getLiftingSurfs()[i]->getWakePanels();
        for (int j=0; j<tempW.size(); j++)
        {
            wPanels.push_back(tempW[j]);
        }
    }
    
    unsigned long nBodyPans = bPanels.size();
    unsigned long nWakePans = wPanels.size();
    
    Ab.resize(nBodyPans,nBodyPans);
    double influence;
    double interpCoeff;
    Eigen::MatrixXi indices(nWakePans,4);
    doubletInfs.resize(nBodyPans);
    sourceInfs.resize(nBodyPans);
    firstDists.resize(nBodyPans);

    std::cout << "Computing Doublet Influence Coefficients..." << std::endl;
    
    for (int i=0; i<nBodyPans; i++)
    {
        for (int j=0; j<nBodyPans; j++)
        {
            if (i==j)
            {
                Ab(i,j) = -0.5;
            }
            else
            {
//                std::cout << j << ", ";
                Ab(i,j) = bPanels[j]->doubletPhi(1,bPanels[i]->getCenter()-0.000001*bPanels[i]->getNormal());
            }
            if (i==0)
            {
                doubletInfs(j) = Ab(i,j);
                firstDists(j) = (bPanels[i]->getCenter()-bPanels[j]->getCenter()).norm();
            }
            
//            std::cout << Ab(i,j) << std::endl;
        }
        for (int j=0; j<nWakePans; j++)
        {
            
            std::vector<bodyPanel*> interpPans; // [Upper1 Lower1 Upper2 Lower2]  Panels that start the bounding wakelines of the wake panel.  Doublet strength is constant along wakelines (muUpper-muLower) and so the doublet strength used for influence of wake panel is interpolated between wakelines.
            
            influence = wPanels[j]->influenceCoeffPhi(bPanels[i]->getCenter(),interpPans,interpCoeff);
            
            
            if (i==0)
            {
                indices.row(j) = getIndices(interpPans);
            }
            Ab(i,indices(j,0)) += influence*(1-interpCoeff);
            Ab(i,indices(j,1)) += influence*(interpCoeff-1);
            Ab(i,indices(j,2)) += influence*interpCoeff;
            Ab(i,indices(j,3)) -= influence*interpCoeff;
        }
    }
    Eigen::VectorXd RHS = getRHS();
    Eigen::VectorXd doubletStrengths(bPanels.size());
    Eigen::VectorXd sourceStrengths(bPanels.size());
    
    doubletStrengths = Ab.householderQr().solve(RHS);
    
    Eigen::VectorXd potentials(nBodyPans);
    Eigen::MatrixXd centers(nBodyPans,3);
    for (int i=0; i<nBodyPans; i++)
    {
        bPanels[i]->setMu(doubletStrengths(i));
        bPanels[i]->setPotential(Vinf);
        potentials(i) = bPanels[i]->getPotential();
//        if (bPanels[i]->getCenter()(0) < (-.95))
//        {
//            std::cout << bPanels[i]->getCenter()(0) << "," << doubletStrengths(i) << "," << potentials(i) << std::endl;
//        }
        sourceStrengths(i) = bPanels[i]->getSigma();
        centers.row(i) = bPanels[i]->getCenter();
    }
    std::cout << doubletStrengths.maxCoeff() << std::endl;
    std::cout << doubletStrengths.minCoeff() << std::endl;
    std::cout << sourceStrengths.maxCoeff() << std::endl;
    std::cout << sourceStrengths.minCoeff() << std::endl;
    
    
    for (int i=0; i<nWakePans; i++)
    {
        wPanels[i]->setMu();
    }
    
    writeVTU(outFile);
}

Eigen::Vector4i runCase::getIndices(std::vector<bodyPanel*> interpPans)
{
    Eigen::Vector4i indices;
    for (int i=0; i<interpPans.size(); i++)
    {
        for (int j=0; j<bPanels.size(); j++)
        {
            if (interpPans[i] == bPanels[j])
            {
                indices(i) = j;
                break;
            }
        }
    }
    return indices;
}

Eigen::VectorXd runCase::getRHS()
{
    sourceInfs(0) = 0;
    Eigen::VectorXd RHS = Eigen::VectorXd::Zero(bPanels.size());;
    Bb.resize(bPanels.size(),bPanels.size());
    for (int i=0; i<bPanels.size(); i++)
    {
        for (int j=0; j<bPanels.size(); j++)
        {
            if (i==j)
            {
//                RHS(i) += 0.5*bPanels[j]->getSigma();
//                RHS(i) += bPanels[j]->sourcePhi(bPanels[j]->getSigma(), bPanels[i]->getCenter());
                RHS(i) += 0;
            }
            else
            {
                if (i==0)
                {
//                    std::cout << j << ", ";
                    sourceInfs(j) = bPanels[j]->sourcePhi(1, bPanels[i]->getCenter());
                    RHS(i) += bPanels[j]->getSigma()*sourceInfs(j);

                }
                else
                {
                    RHS(i) += bPanels[j]->sourcePhi(bPanels[j]->getSigma(), bPanels[i]->getCenter());
                }
            }
//            std::cout << RHS(i) << std::endl;
        }
    }
    return -RHS;
}

void runCase::writeVTU(std::string filename)
{
    std::ofstream fid;
    fid.open(filename);
    if (fid.is_open())
    {
        fid << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << std::endl;
        fid << "  <UnstructuredGrid>" << std::endl;
        fid << "    <Piece NumberOfPoints=\"" << geom->getNumberOfNodes() << "\" NumberOfCells=\"" << bPanels.size() << "\">" << std::endl;
        fid << "      <CellData Scalars=\"scalars\">" << std::endl;
        fid << "        <DataArray type=\"Float64\" Name=\"Potential\" NumberOfComponents=\"1\" Format=\"ascii\">" << std::endl;
        for (int i=0; i<bPanels.size(); i++)
        {
            fid << bPanels[i]->getPotential() << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "        <DataArray type=\"Float64\" Name=\"DoubletStrength\" NumberOfComponents=\"1\" Format=\"ascii\">" << std::endl;
        for (int i=0; i<bPanels.size(); i++)
        {
            fid << bPanels[i]->getMu() << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "        <DataArray type=\"Float64\" Name=\"DoubletInfluence\" NumberOfComponents=\"1\" Format=\"ascii\">" << std::endl;
        for (int i=0; i<bPanels.size(); i++)
        {
            fid << doubletInfs(i) << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "        <DataArray type=\"Float64\" Name=\"SourceInfluence\" NumberOfComponents=\"1\" Format=\"ascii\">" << std::endl;
        for (int i=0; i<bPanels.size(); i++)
        {
            fid << sourceInfs(i) << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "        <DataArray type=\"Float64\" Name=\"Distance\" NumberOfComponents=\"1\" Format=\"ascii\">" << std::endl;
        for (int i=0; i<bPanels.size(); i++)
        {
            fid << firstDists(i) << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "      </CellData>" << std::endl;
        fid << "      <Points>" << std::endl;
        fid << "        <DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
        Eigen::MatrixXd nodes = geom->getNodes();
        for (int i=0; i<nodes.rows(); i++)
        {
            fid << nodes(i,0) << "  " << nodes(i,1) << "  " << nodes(i,2) << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "      </Points>" << std::endl;
        fid << "      <Cells>" << std::endl;
        fid << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">" << std::endl;
        Eigen::VectorXi verts;
        for (int i=0; i<bPanels.size(); i++)
        {
            verts = bPanels[i]->getVerts();
            fid << verts(0) << "  " << verts(1) << "  " << verts(2) << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">" << std::endl;
        for (int i=0; i<bPanels.size(); i++)
        {
            fid << 3*(i+1) << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "        <DataArray type=\"UInt8\" Name=\"types\" Format=\"ascii\">" << std::endl;
        for (int i=0; i<bPanels.size(); i++)
        {
            fid << 5 << std::endl;
        }
        fid << "        </DataArray>" << std::endl;
        fid << "      </Cells>" << std::endl;
        fid << "    </Piece>" << std::endl;
        fid << "  </UnstructuredGrid>" << std::endl;
        fid << "</VTKFile>" << std::endl;
        std::cout << "Data written to " << filename << std::endl;
    }
    fid.close();
}
