//
//  runCase.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "cpCase.h"

cpCase::~cpCase()
{
    for (int i=0; i<bStreamlines.size(); i++)
    {
        delete bStreamlines[i];
    }
}

void cpCase::run()
{
    bool converged;
    std::string check = "\u2713";
    setSourceStrengths();
    converged = solveMatrixEq();
    std::cout << std::setw(17) << std::left << check << std::flush;
    compVelocity();
    std::cout << std::setw(17) << std::left << check << std::flush;
    trefftzPlaneAnalysis();
    std::cout << std::setw(18) << std::left << check << std::flush;

    createStreamlines();
    std::cout << std::setw(16) << std::left << check << std::endl;
        
    if (!converged)
    {
        std::cout << "*** Warning : Solution did not converge ***" << std::endl;;
    }
    
    writeFiles();
}

Eigen::Vector3d cpCase::windToBody(double V, double alpha, double beta)
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
    transform = T;

    Vel = transform*Vt;
    
    return Vel;
}

Eigen::Vector3d cpCase::bodyToWind(const Eigen::Vector3d &vec)
{
    Eigen::Matrix3d T = transform.transpose();
    return T*vec;
}

void cpCase::setSourceStrengths()
{
    sigmas.resize(bPanels->size());
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setSigma(Vinf,0);
        sigmas(i) = (*bPanels)[i]->getSigma();
    }
}

bool cpCase::solveMatrixEq()
{
    bool converged = true;
    
    // Solve matrix equations and set potential for all panels;
    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::VectorXd RHS = -(*B)*sigmas;
    Eigen::VectorXd doubletStrengths(bPanels->size());
  
    
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute((*A));
    doubletStrengths = res.solve(RHS);
    if (res.error() > pow(10,-10))
    {
        converged = false;
    }
    
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(i));
        (*bPanels)[i]->setPotential(Vinf);
    }
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf);
    }
    return converged;
}

void cpCase::compVelocity()
{
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;
    for (int i=0; i<bPanels->size(); i++)
    {
        p = (*bPanels)[i];
        p->computeVelocity();
        p->computeCp(Vmag,PG);
        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
        moment = p->computeMoments(params->cg);
        CM(0) += moment(0)/(params->Sref*params->bref);
        CM(1) += moment(1)/(params->Sref*params->cref);
        CM(2) += moment(2)/(params->Sref*params->bref);
    }
    Fwind = bodyToWind(Fbody);
}

void cpCase::trefftzPlaneAnalysis()
{
    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    for (int i=0; i<wakes.size(); i++)
    {
        wakes[i]->trefftzPlane(Vmag,params->Sref);
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
}

void cpCase::createStreamlines()
{
    // Gather TE Panels
    std::vector<surface*> surfs = geom->getSurfaces();
    std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> streamPnts;
    bodyStreamline* s;
    
    for (int i=0; i<surfs.size(); i++)
    {
        streamPnts = surfs[i]->getStreamlineStartPnts(Vinf);
        for (int j=0; j<streamPnts.size(); j++)
        {
            s = new bodyStreamline(std::get<0>(streamPnts[j]),std::get<1>(streamPnts[j]),Vinf,geom,2,false);
            bStreamlines.push_back(s);
        }
    }
}

void cpCase::writeFiles()
{
    std::stringstream caseLabel;
    caseLabel << "/V" << Vmag << "_Mach" << mach << "_alpha" << alpha << "_beta" << beta;
    boost::filesystem::path subdir = boost::filesystem::current_path().string()+caseLabel.str();
    if (!boost::filesystem::exists(subdir))
    {
        boost::filesystem::create_directories(subdir);
    }
    Eigen::MatrixXd nodeMat = geom->getNodePnts();
    writeBodyData(subdir,nodeMat);
    if (geom->getWakes().size() > 0)
    {
        writeWakeData(subdir,nodeMat);
        writeSpanwiseData(subdir);
    }
    
    writeBodyStreamlines(subdir);
}

void cpCase::writeBodyData(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals");
    Eigen::MatrixXi con(bPanels->size(),3);
    mu.data.resize(bPanels->size(),1);
    pot.data.resize(bPanels->size(),1);
    V.data.resize(bPanels->size(),3);
    Cp.data.resize(bPanels->size(),1);
    bN.data.resize(bPanels->size(),3);
    for (int i=0; i<bPanels->size(); i++)
    {
        mu.data(i,0) = (*bPanels)[i]->getMu();
        pot.data(i,0) = (*bPanels)[i]->getPotential();
        V.data.row(i) = (*bPanels)[i]->getGlobalV();
        Cp.data(i,0) = (*bPanels)[i]->getCp();
        con.row(i) = (*bPanels)[i]->getVerts();
        bN.data.row(i) = (*bPanels)[i]->getBezNormal();
    }
    
    data.push_back(mu);
    data.push_back(pot);
    data.push_back(V);
    data.push_back(Cp);
    data.push_back(bN);
    
    piece body;
    body.pnts = nodeMat;
    body.connectivity = con;
    body.cellData = data;
    
    std::string fname = path.string()+"/surfaceData.vtu";
    VTUfile bodyFile(fname,body);
}

void cpCase::writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::MatrixXi con(wPanels->size(),3);
    mu.data.resize(wPanels->size(),1);
    pot.data.resize(wPanels->size(),1);
    for (int i=0; i<wPanels->size(); i++)
    {
        mu.data(i,0) = (*wPanels)[i]->getMu();
        pot.data(i,0) = (*wPanels)[i]->getPotential();
        con.row(i) = (*wPanels)[i]->getVerts();
    }
    
    data.push_back(mu);
    data.push_back(pot);
    
    piece wake;
    wake.pnts = nodeMat;
    wake.connectivity = con;
    wake.cellData = data;
    
    std::string fname = path.string()+"/wakeData.vtu";
    VTUfile wakeFile(fname,wake);
}

void cpCase::writeSpanwiseData(boost::filesystem::path path)
{
    std::vector<wake*> wakes = geom->getWakes();
    for (int i=0; i<wakes.size(); i++)
    {
        Eigen::VectorXd spanLoc,Cl,Cd;
        spanLoc = 2*wakes[i]->getSpanwisePnts()/params->Sref;
        Cl = wakes[i]->getSpanwiseCl()/PG;
        Cd = wakes[i]->getSpanwiseCd()/pow(PG,2);
        
        std::stringstream ss;
        ss << path.string() << "/spanwiseData_Wake" << i+1 << ".csv";
        std::string fname = ss.str();
        std::ofstream fout;
        fout.open(fname);
        if (fout)
        {
            fout << "2y/b,Cl,Cdi" << std::endl;
            for (int i=0; i<spanLoc.size(); i++)
            {
                fout << spanLoc(i) << "," << Cl(i) << "," << Cd(i) << std::endl;
            }
        }
        fout.close();
    }
}

void cpCase::writeBodyStreamlines(boost::filesystem::path path)
{
    piece p;
    std::vector<piece> pieces;
    std::vector<pntDataArray> data;
    pntDataArray vel("Velocity");
    Eigen::MatrixXi con;
    Eigen::MatrixXd pntMat;
    std::vector<Eigen::Vector3d> pnts,velocities;
    
    for (int i=0; i<bStreamlines.size(); i++)
    {
        pnts = bStreamlines[i]->getPnts();
        velocities = bStreamlines[i]->getVelocities();
        vel.data.resize(velocities.size(),3);
        pntMat.resize(pnts.size(),3);
        con.resize(pnts.size()-1,2);
        for (int j=0; j<pnts.size(); j++)
        {
            pntMat.row(j) = pnts[j];
            vel.data.row(j) = velocities[j];
            if (j<con.rows())
            {
                con(j,0) = j;
                con(j,1) = j+1;
            }
        }
        data.push_back(vel);
        p.pnts = pntMat;
        p.connectivity = con;
        p.pntData = data;
        
        pieces.push_back(p);
        data.clear();
    }
    
    
    
    std::string fname = path.string()+"/streamlines.vtu";
    VTUfile wakeFile(fname,pieces);
}


