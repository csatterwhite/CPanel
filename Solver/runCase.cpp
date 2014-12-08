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
        }
    }
}

void runCase::solveMatrixEq()
{
    
    // Get panels into vectors of body panels and wake panels. Trailing edge panels are placed at the beginning of the body panel vector to speed up search for the body panels shedding wake lines in the wake influence calculation.
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
    
    // Construct doublet and source influence coefficient matrices for body panels
    unsigned long nBodyPans = bPanels.size();
    unsigned long nWakePans = wPanels.size();
    
    Ab.resize(nBodyPans,nBodyPans);
    Bb.resize(nBodyPans,nBodyPans);
    sigmas.resize(nBodyPans);

    std::cout << "Computing Influence Coefficients..." << std::endl;
    Eigen::Vector3i percentage;
    percentage << 25,50,75;
    
    for (int j=0; j<nBodyPans; j++)
    {
        for (int i=0; i<nBodyPans; i++)
        {
            bPanels[j]->panelPhiInf(bPanels[i]->getCenter(),Bb(i,j),Ab(i,j));
            if (i==j)
            {
                Ab(i,j) = -0.5;
            }
            if (j==0)
            {
                sigmas(i) = bPanels[i]->getSigma();
            }
        }
        for (int i=0; i<percentage.size(); i++)
        {
            if ((100*j/nBodyPans) <= percentage(i) && 100*(j+1)/nBodyPans > percentage(i))
            {
                std::cout << percentage(i) << "%" << std::endl;
            }
        }
    }
    std::cout << "Complete" << std::endl;
    // Add influence of wake panels to doublet influence coefficient matrix
    
    std::vector<bodyPanel*> interpPans(4); // [Upper1 Lower1 Upper2 Lower2]  Panels that start the bounding wakelines of the wake panel.  Doublet strength is constant along wakelines (muUpper-muLower) and so the doublet strength used for influence of wake panel is interpolated between wakelines.
    double interpCoeff;
    double influence;
    Eigen::Vector4i indices;
    
    std::cout << "Computing wake panel influences..." << std::endl;
    for (int j=0; j<nWakePans; j++)
    {
        wPanels[j]->interpPanels(interpPans,interpCoeff);
        indices = getIndices(interpPans);
        for (int i=0; i<nBodyPans; i++)
        {
            influence = wPanels[j]->dubPhiInf(bPanels[i]->getCenter());
            Ab(i,indices(0)) += influence*(1-interpCoeff);
            Ab(i,indices(1)) += influence*(interpCoeff-1);
            Ab(i,indices(2)) += influence*interpCoeff;
            Ab(i,indices(3)) -= influence*interpCoeff;
        }
        for (int i=0; i<percentage.size(); i++)
        {
            if ((100*j/nWakePans) <= percentage(i) && 100*(j+1)/nWakePans > percentage(i))
            {
                std::cout << percentage(i) << "%" << std::endl;
            }
        }
        
    }
    std::cout << "Complete" << std::endl;
    
    // Solve matrix equations and set potential for all panels;
    
    Eigen::VectorXd RHS = -Bb*sigmas;
    Eigen::VectorXd doubletStrengths(bPanels.size());
    
    std::cout << "Solving system of equations..." << std::endl;
    time_t ts;
    time(&ts);
    time_t tf;
    
//    Eigen::ConjugateGradient<Eigen::MatrixXd> cg(Ab);
//    doubletStrengths = cg.solve(RHS);
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute(Ab);
    doubletStrengths = res.solve(RHS);
    std::cout << "#iterations:     " << res.iterations() << std::endl;
    std::cout << "estimated error: " << res.error()      << std::endl;
//    doubletStrengths = Ab.householderQr().solve(RHS);
    time(&tf);
    std::cout << "Time to solve Ax=b : " << difftime(tf,ts) << "seconds" << std::endl;
    
    std::cout << "Complete" << std::endl;
    for (int i=0; i<nBodyPans; i++)
    {
        bPanels[i]->setMu(doubletStrengths(i));
        bPanels[i]->setPotential(Vinf);
    }
    for (int i=0; i<nWakePans; i++)
    {
        wPanels[i]->setMu();
        wPanels[i]->setPotential(Vinf);
    }
    
    //  Velocity Survey with known doublet and source strengths
    
    std::cout << "Calculating Velocities..." << std::endl;
    allPanels.insert(allPanels.end(),bPanels.begin(),bPanels.end());
    allPanels.insert(allPanels.end(),wPanels.begin(),wPanels.end());
//    Eigen::MatrixXd velocities = Eigen::MatrixXd::Zero(allPanels.size(),3);
//
//    
//    for (int i=0; i<allPanels.size(); i++)
//    {
//        allPanels[i]->computeVelocity();
//        velocities.row(i) = allPanels[i]->getGlobalV();
//        std::cout << i << std::endl;
//        for (int j=0; j<percentage.size(); j++)
//        {
//            if ((100*i/allPanels.size()) <= percentage(j) && 100*(i+1)/allPanels.size() > percentage(j))
//            {
//                std::cout << percentage(j) << "%" << std::endl;
//            }
//        }
//    }
    
    Eigen::MatrixXd velocities = Eigen::MatrixXd::Zero(bPanels.size(),3);
    
    
    for (int i=0; i<bPanels.size(); i++)
    {
        bPanels[i]->computeVelocity();
        velocities.row(i) = bPanels[i]->getGlobalV();
        for (int j=0; j<percentage.size(); j++)
        {
            if ((100*i/bPanels.size()) <= percentage(j) && 100*(i+1)/bPanels.size() > percentage(j))
            {
                std::cout << percentage(j) << "%" << std::endl;
            }
        }
    }
    
    std::cout << "Complete" << std::endl;

    std::cout << "Writing .vtu file" << std::endl;
    
    std::vector<cellDataArray*> bCellData,wCellData;
    cellDataArray bMu,bPot,bV,bNeighbs,wMu,wPot,wV,wNeighbs;
    Eigen::MatrixXd bDub(bPanels.size(),1),bPotential(bPanels.size(),1),bVel(bPanels.size(),3),bNeighbors(bPanels.size(),1),wDub(wPanels.size(),1),wPotential(wPanels.size(),1),wVel(wPanels.size(),3),wNeighbors(wPanels.size(),1);
    Eigen::MatrixXi bodyCon(bPanels.size(),3),wakeCon(wPanels.size(),3);
    for (int i=0; i<bPanels.size(); i++)
    {
        bDub(i,0) = bPanels[i]->getMu();
        bPotential(i,0) = bPanels[i]->getPotential();
        bVel.row(i) = bPanels[i]->getGlobalV();
        bodyCon.row(i) = bPanels[i]->getVerts();
        bNeighbors(i,0) = bPanels[i]->getNeighbors().size();
    }
    for (int i=0; i<wPanels.size(); i++)
    {
        wDub(i,0) = wPanels[i]->getMu();
        wPotential(i,0) = wPanels[i]->getPotential();
//        wVel.row(i) = wPanels[i]->getGlobalV();
        wakeCon.row(i) = wPanels[i]->getVerts();
        wNeighbors(i,0) = wPanels[i]->getNeighbors().size();
    }
    bMu.name = "DoubletStrengths";
    bMu.data = bDub;
    bCellData.push_back(&bMu);
    
    bPot.name = "Velocity Potential";
    bPot.data = bPotential;
    bCellData.push_back(&bPot);
    
//    bV.name = "Velocity";
//    bV.data = bVel;
//    bCellData.push_back(&bV);
    
    bNeighbs.name = "Number of Neighbors";
    bNeighbs.data = bNeighbors;
    bCellData.push_back(&bNeighbs);
    
    piece body;
    body.pnts = geom->getNodes();
    body.connectivity = bodyCon;
    body.cellData = bCellData;
    
    wMu.name = "DoubletStrengths";
    wMu.data = wDub;
    wCellData.push_back(&wMu);
    
    wPot.name = "Velocity Potential";
    wPot.data = wPotential;
    wCellData.push_back(&wPot);
    
    wNeighbs.name = "Number of Neighbors";
    wNeighbs.data = wNeighbors;
    wCellData.push_back(&wNeighbs);
    
    piece wake;
    wake.pnts = geom->getNodes();
    wake.connectivity = wakeCon;
    wake.cellData = wCellData;
    
    std::vector<piece*> pieces;
    pieces.push_back(&body);
    pieces.push_back(&wake);
    VTUfile VTU(pieces);
    VTU.write(outFile);
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
