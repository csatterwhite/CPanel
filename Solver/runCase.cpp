//
//  runCase.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "runCase.h"

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
            if (tempB[j]->isUpper() || tempB[j]->isLower())
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
    
////    std::vector<horseshoeVortex*> hvec;
//    std::vector<wakePanel*> vs;
//    std::vector<wake*> wakes = geom->getWakes();
//    for (int i=0; i<wakes.size(); i++)
//    {
////        std::vector<horseshoeVortex*> hnew = wakes[i]->getHorseshoes();
////        hvec.insert(hvec.end(), hnew.begin(), hnew.end());
//        std::vector<wakePanel*> vnew = wakes[i]->getVortexSheets();
//        vs.insert(vs.end(),vnew.begin(),vnew.end());
//    }
//    bodyPanel* upper;
//    bodyPanel* lower;
//    double influence;
//    int uIndex,lIndex;
////    for (int i=0; i<hvec.size(); i++)
//    for (int i=0; i<vs.size(); i++)
//    {
////        upper = hvec[i]->getUpper();
////        lower = hvec[i]->getLower();
//        upper = vs[i]->getUpper();
//        lower = vs[i]->getLower();
//        // Find indices of upper and lower panels
//        uIndex = -1;
//        lIndex = -1;
//        for (int j=0; j<bPanels.size(); j++)
//        {
//            if (bPanels[j]==upper)
//            {
//                uIndex = j;
//            }
//            else if (bPanels[j]==lower)
//            {
//                lIndex = j;
//            }
//            if (uIndex > 0 && lIndex > 0)
//            {
//                break;
//            }
//        }
//        for (int j=0; j<bPanels.size(); j++)
//        {
////            influence = hvec[i]->phiInfluence(bPanels[j]->getCenter());
//            influence = vs[i]->dubPhiInf(bPanels[j]->getCenter());
//            Ab(j,uIndex) += influence;
//            Ab(j,lIndex) -= influence;
//        }
//    }
    
    
    
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
    
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute(Ab);
    doubletStrengths = res.solve(RHS);
    std::cout << "#iterations:     " << res.iterations() << std::endl;
    std::cout << "estimated error: " << res.error()      << std::endl;
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
//    for (int i=0; i<hvec.size(); i++)
//    {
//        hvec[i]->setStrength();
//    }
//    for (int i=0; i<vs.size(); i++)
//    {
//        vs[i]->setStrength();
//    }
    
    //  Velocity Survey with known doublet and source strengths
    Eigen::Vector3d BodyF = Eigen::Vector3d::Zero();
    std::cout << "Calculating Velocities..." << std::endl;
    for (int i=0; i<bPanels.size(); i++)
    {
        bPanels[i]->computeVelocity();
        bPanels[i]->computeCp(Vmag);
        BodyF += -bPanels[i]->getCp()*bPanels[i]->getArea()*bPanels[i]->getNormal()/6; // 6 is Sref
        for (int j=0; j<percentage.size(); j++)
        {
            if ((100*i/bPanels.size()) <= percentage(j) && 100*(i+1)/bPanels.size() > percentage(j))
            {
                std::cout << percentage(j) << "%" << std::endl;
            }
        }
    }
    std::cout << "Complete" << std::endl;
//    std::cout << "Body Forces\t" << BodyF(0) << "\t" << BodyF(1) << "\t" << BodyF(2) << std::endl;
    
    /////////// TESTING POTENTIAL FROM WAKE /////////
    
//    int n = 20;
//    Eigen::MatrixXd testPnts(n,3);
//    Eigen::VectorXd testPhi = Eigen::VectorXd::Zero(n);
//    for (int i=0; i<n; i++)
//    {
////        testPnts(i,0) = 20+i*(1.0/(n-1));
//        testPnts(i,0) = 20;
//        testPnts(i,1) = 2.72;
//        testPnts(i,2) = 0.001+i*(0.06/(n-1));
////        testPnts(i,2) = -0.1;
//        std::vector<wakePanel*> vs = geom->getWakes()[0]->getVortexSheets();
//        for (int j=0; j<vs.size(); j++)
//        {
//            if (i==0)
//            {
//                std::cout << vs[j]->getCenter()(1) << "\t" << vs[j]->getMu() << std::endl;
//            }
//            testPhi(i) += vs[j]->panelPhi(testPnts.row(i));
//        }
////        if (testPnts(i,2) > 0)
////        {
////            testPhi(i) -= wStrength;
////        }
//        std::cout << testPnts(i,0) << "\t" << testPnts(i,1) << "\t" << testPnts(i,2) << "\t" << testPhi(i) << std::endl;
//    }
//    
    
    /////////////////////////////////////////////////
    
    double CL,CD;
    Eigen::VectorXd yLoc,Cl,Cd;
    Eigen::MatrixXd Fsect;
    Eigen::Vector3d Fbody;
    std::vector<wake*> wakes = geom->getWakes();
    for (int i=0; i<wakes.size(); i++)
    {
//        wakes[i]->trailingEdge(Vinf, 6, Fbody, yLoc, Fsect, 1);
        wakes[i]->trefftzPlane(Vmag,6,CL,CD,yLoc,Cl,Cd);
//        wakes[i]->horseshoeTrefftz(Vmag, 6, CL, CD, yLoc, Cl, Cd);
//        wakes[i]->sheetTrefftz(Vinf.norm(), 6, CL, CD, yLoc, Cl, Cd);
    }
//    L = Fbody(2)*cos(alpha*M_PI/180)-Fbody(0)*sin(alpha*M_PI/180);
//    D = Fbody(2)*sin(alpha*M_PI/180)+Fbody(0)*cos(alpha*M_PI/180);
    std::cout << "CL = " << CL << "\t\tCD = " << CD << std::endl;
    for (int i=0; i<yLoc.rows(); i++)
    {
        std::cout << yLoc(i) << "\t" << Cl(i) << "\t" << Cd(i) << std::endl;
    }

    std::cout << "Writing .vtu file" << std::endl;
    
    writeFiles();
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

void runCase::writeFiles()
{
    std::stringstream caseLabel;
    caseLabel << "/V" << Vmag << "_alpha" << alpha << "_beta" << beta;
    boost::filesystem::path subdir = path+geomName+caseLabel.str();
    if (!boost::filesystem::exists(subdir))
    {
        boost::filesystem::create_directories(subdir);
    }
    writeBodyData(subdir);
//    if (wPanels.size() > 0)
//    {
//        writeWakeData(subdir);
//    }
}

void runCase::writeBodyData(boost::filesystem::path path)
{
    std::vector<cellDataArray*> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp");
    Eigen::MatrixXi con(bPanels.size(),3);
    mu.data.resize(bPanels.size(),1);
    pot.data.resize(bPanels.size(),1);
    V.data.resize(bPanels.size(),3);
    Cp.data.resize(bPanels.size(),1);
    for (int i=0; i<bPanels.size(); i++)
    {
        mu.data(i,0) = bPanels[i]->getMu();
        pot.data(i,0) = bPanels[i]->getPotential();
        V.data.row(i) = bPanels[i]->getGlobalV();
        Cp.data(i,0) = bPanels[i]->getCp();
        con.row(i) = bPanels[i]->getVerts();
    }
    
    data.push_back(&mu);
    data.push_back(&pot);
    data.push_back(&V);
    data.push_back(&Cp);
    
    piece body;
    body.pnts = geom->getNodes();
    body.connectivity = con;
    body.cellData = data;
    
    std::string fname = path.string()+"/body.vtu";
    VTUfile bodyFile(fname,&body);
}

void runCase::writeWakeData(boost::filesystem::path path)
{
    std::vector<cellDataArray*> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::MatrixXi con(wPanels.size(),3);
    mu.data.resize(wPanels.size(),1);
    pot.data.resize(wPanels.size(),1);
    for (int i=0; i<wPanels.size(); i++)
    {
        mu.data(i,0) = wPanels[i]->getMu();
        pot.data(i,0) = wPanels[i]->getPotential();
        con.row(i) = wPanels[i]->getVerts();
    }
    
    data.push_back(&mu);
    data.push_back(&pot);
    
    piece wake;
    wake.pnts = geom->getNodes();
    wake.connectivity = con;
    wake.cellData = data;
    
    std::string fname = path.string()+"/wake.vtu";
    VTUfile wakeFile(fname,&wake);
}
