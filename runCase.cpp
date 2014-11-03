//
//  runCase.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "runCase.h"

runCase::runCase(geometry *geom,double V,double alpha,double beta) : geom(geom)
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
    std::vector<surface*> surfaces = geom->getSurfaces();
    std::vector<bodyPanel*> bPanels;
    std::vector<wakePanel*> wPanels;
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

    std::cout << "Computing Doublet Influence Coefficients..." << std::endl;
    Eigen::VectorXi percentComplete(4);
    percentComplete(0) = 25;
    percentComplete(1) = 50;
    percentComplete(2) = 75;
    percentComplete(3) = 100;
    
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
                Ab(i,j) = bPanels[j]->doubletPhi(1,bPanels[i]->getCenter());
            }
        }
        for (int j=0; j<nWakePans; j++)
        {
            
            std::vector<bodyPanel*> interpPans; // [Upper1 Lower1 Upper2 Lower2]  Panels that start the bounding wakelines of the wake panel.  Doublet strength is constant along wakelines (muUpper-muLower) and so the doublet strength used for influence of wake panel is interpolated between wakelines.
            
            influence = wPanels[j]->influenceCoeffPhi(bPanels[i]->getCenter(),interpPans,interpCoeff);
            
            
            if (i==0)
            {
                indices.row(j) = getIndices(interpPans,bPanels);
            }
            Ab(i,indices(j,0)) += influence*(1-interpCoeff);
            Ab(i,indices(j,1)) += influence*(interpCoeff-1);
            Ab(i,indices(j,2)) += influence*interpCoeff;
            Ab(i,indices(j,3)) -= influence*interpCoeff;
        }
        if (i>0)
        {
            for (int k=0; k<percentComplete.size(); k++)
            {
                if (i/nBodyPans >= (percentComplete(k)/100) && (i-1)/nBodyPans < (percentComplete(k)/100))
                {
                    std::cout << percentComplete(k) << "%" << std::endl;
                }
            }
        }
    }
    std::cout << std::endl;
    Eigen::VectorXd RHS = getRHS(bPanels);
    Eigen::VectorXd doubletStrengths(bPanels.size());
    Eigen::VectorXd sourceStrengths(bPanels.size());
    
    doubletStrengths = Ab.householderQr().solve(RHS);
    
    Eigen::VectorXd potentials(nBodyPans);
    Eigen::MatrixXd centers(nBodyPans,3);
    for (int i=0; i<nBodyPans; i++)
    {
        bPanels[i]->setMu(doubletStrengths(i));
        bPanels[i]->setPotential(Vinf);
        potentials[i] = bPanels[i]->getPotential();
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
    
}

Eigen::Vector4i runCase::getIndices(std::vector<bodyPanel*> interpPans, std::vector<bodyPanel*> bPanels)
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

Eigen::VectorXd runCase::getRHS(const std::vector<bodyPanel*> &bPanels)
{
    Eigen::VectorXd RHS(bPanels.size());
    Bb.resize(bPanels.size(),bPanels.size());
    for (int i=0; i<bPanels.size(); i++)
    {
        RHS(i) = bPanels[i]->getSigma();
        for (int j=0; j<bPanels.size(); j++)
        {
            if (i==j)
            {
                Bb(i,j) = 0;
            }
            else
            {
                Bb(i,j) = bPanels[j]->sourcePhi(1,bPanels[i]->getCenter());
            }
        }
    }
    RHS = -Bb*RHS;
    return RHS;
}
