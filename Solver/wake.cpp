//
//  wake.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wake.h"

wake::~wake()
{
    for (int i=0; i<wpanels.size(); i++)
    {
        delete wpanels[i];
    }
    for (int i=0; i<wakeLines.size(); i++)
    {
        delete wakeLines[i];
    }
    for (int i=0; i<horseShoes.size(); i++)
    {
        delete horseShoes[i];
    }
    for (int i=0; i<vortexSheets.size(); i++)
    {
        delete vortexSheets[i];
    }
}

wake::wake(const wake& copy)
{
    for (int i=0; i<copy.wpanels.size(); i++)
    {
        wpanels[i] = new wakePanel(*copy.wpanels[i]);
    }
}

void wake::addPanel(wakePanel* wPan)
{
    Eigen::Vector3i verts = wPan->getVerts();
    Eigen::Vector3d pnt;
    if (wpanels.size() == 0)
    {
        // Initialize wake dimensions with first panel
        yMax = nodes->row(verts(0))(1);
        yMin = yMax;
        x0 = nodes->row(verts(0))(0);
        xf = x0;
    }
    
    for (int i=0; i<verts.size(); i++)
    {
        pnt = nodes->row(verts(i));
        if (pnt(1) > yMax)
        {
            yMax = pnt(1);
        }
        else if (pnt(1) < yMin)
        {
            yMin = pnt(1);
        }
        if (pnt(0) > xf)
        {
            xf = pnt(0);
            zf = pnt(2);
        }
        else if (pnt(0) < x0)
        {
            x0 = pnt(0);
            z0 = pnt(2);
        }
    }
    wpanels.push_back(wPan);
}

void wake::setNeighbors(panelOctree* oct)
{
    for (int i=0; i<wpanels.size(); i++)
    {
        short nVerts = wpanels[i]->getVerts().size();
        if (edgeVerts(wpanels[i]) == 2)
        {
            // Panel is on edge of wake and will have two neighbors for tris and three for quads
            wpanels[i]->setNeighbors(oct,nVerts-1);
        }
        else if (edgeVerts(wpanels[i]) == 3)
        {
            // Panel is in downstream corner and will have 1 neighbor for tris and 2 for quads
            wpanels[i]->setNeighbors(oct,nVerts-2);
        }
        else
        {
            wpanels[i]->setNeighbors(oct,nVerts);
        }
        if (wpanels[i]->isTEpanel())
        {
            wpanels[i]->setParentPanels();
            wpanels[i]->addWakeLine(wakeLines);
            horseshoeVortex* h = wpanels[i]->makeHorseshoe();
            horseShoes.push_back(h);
            wakePanel* w = wpanels[i]->makeVortexSheet();
            vortexSheets.push_back(w);
        }
    }
}

short wake::edgeVerts(wakePanel* p)
{
    // Returns number of verts on edge of wake. 2 indicates panel edge of wake. 3 indicates panel in corner at downstream side of wake;
    short count = 0;
    Eigen::VectorXi verts = p->getVerts();
    Eigen::Vector3d pnt;
    for (int i=0; i<verts.size(); i++)
    {
        pnt = nodes->row(verts(i));
        if (pnt(0) == xf || pnt(1) == yMin || pnt(1) == yMax)
        {
            count++;
        }
    }
    
    return count;
}

void wake::trefftzPlane(double Vinf,double Sref, double &CL, double &CD, Eigen::VectorXd &yLoc, Eigen::VectorXd &Cl, Eigen::VectorXd &Cd)
{
    int nPnts = 10;
    if (nPnts % 2 != 0)
    {
        //Number is odd and needs to be even for simpsons rule integration.
        nPnts++;
    }
    Eigen::VectorXd w,v,dPhi;
    yLoc.resize(nPnts+1);
    yLoc(0) = yMin;
    yLoc(nPnts) = yMax;
    double step = (yMax-yMin)/(nPnts);
    w = Eigen::VectorXd::Zero(nPnts+1);
    dPhi = Eigen::VectorXd::Zero(nPnts+1);
    Cl = Eigen::VectorXd::Zero(nPnts+1);
    Cd = Eigen::VectorXd::Zero(nPnts+1);
    Eigen::MatrixXd trefftzPnts = Eigen::MatrixXd::Zero(nPnts+1,3);

    double xTrefftz = x0+(xf-x0)/2;
    Eigen::Vector3d pWake;
    for (int i=1; i<nPnts; i++)
    {
        yLoc(i) = yMin+i*step;
        pWake = pntInWake(xTrefftz, yLoc(i));
        w(i) = Vradial(pWake);
        dPhi(i) = -wakeStrength(yLoc(i));
        Cl(i) = 2*dPhi(i)/(Vinf*Sref);
        Cd(i) = dPhi(i)*w(i)/(Vinf*Vinf*Sref);
    }
    for (int i=0; i<yLoc.size(); i++)
    {
        std::cout << yLoc(i) << "\t" << w(i) << "\t" << dPhi(i) << "\t" << Cl(i) << "\t" << Cd(i) << std::endl;
    }

    
    int i=0;
    CL = 0;
    CD = 0;
    while (i < Cl.rows()-2)
    {
        CL += 1.0/3*step*(Cl(i)+4*Cl(i+1)+Cl(i+2));
        CD += 1.0/3*step*(Cd(i)+4*Cd(i+1)+Cd(i+2));
        i += 2;
    }
}

void wake::trailingEdge(Eigen::Vector3d Vinf,double Sref, Eigen::Vector3d &Fbody, Eigen::VectorXd &yLoc, Eigen::MatrixXd &Fsect, double rho)
{
    Eigen::VectorXd w,v,dPhi;
    int nLines = (int)wakeLines.size();
    yLoc.resize(nLines);
    Eigen::Vector3d normWake;
    int nSurv = 45;
    Eigen::Matrix3d POIs;
    Eigen::MatrixXd survPnts(nSurv,3);
    Eigen::Matrix<bool,Eigen::Dynamic,1> upperFlag(nSurv,1);
    Eigen::Vector3d vPOI; // Velocity at point just off TE
    Eigen::Vector3d circDir;
    Eigen::Matrix3d Ftemp;
    Fsect.resize(nLines,3);
    double gamma;
    for (int i=0; i<nLines; i++)
    {
        POIs = wakeLines[i]->randSurveyPnts(survPnts,upperFlag,circDir);
        yLoc(i) = wakeLines[i]->getY();
        Eigen::Vector3d vel;
        for (int j=0; j<POIs.rows(); j++)
        {
            vel = pntVel(POIs.row(j), survPnts, upperFlag, Vinf);
            std::cout << POIs(j,1) << "\t" << vel(0) << "\t" << vel(1) << "\t" << vel(2) << std::endl;
            gamma = wakeStrength(POIs(j,1));
            Ftemp.row(j) = rho*(vel.cross(gamma*circDir));
        }
        Fsect.row(i) = 1.0/3*((POIs.row(2)-POIs.row(0)).norm())*(Ftemp.row(0)+4*Ftemp.row(1)+Ftemp.row(2));
        std::cout << "\n" << std::endl;

//        std::cout << yLoc(i) << "\t" << vel(0) << "\t" << vel(1) << "\t" << vel(2) << std::endl;
    }
    
    Fbody = Fsect.colwise().sum();
}


wakeLine* wake::findWakeLine(double y)
{
    double y1,y2;
    for (int i=0; i<wakeLines.size(); i++)
    {
        y1 = wakeLines[i]->getP1()(1);
        y2 = wakeLines[i]->getP2()(1);
        if (y <= y2 && y >= y1)
        {
            return wakeLines[i];
        }
    }
    return nullptr;
}

double wake::wakeStrength(double y)
{
    wakeLine* wl1 = nullptr;
    wakeLine* wl2 = nullptr;
    if (y < wakeLines[1]->getY())
    {
        wl1 = wakeLines[0];
        wl2 = wakeLines[1];
    }
    else if (y >= wakeLines.end()[-1]->getY())
    {
        wl1 = wakeLines.end()[-2]; //Second to last wakeline
        wl2 = wakeLines.end()[-1]; //Last wakeline
    }
    else
    {
        for (int i=1; i<wakeLines.size()-1; i++)
        {
            if ((wakeLines[i]->getY() <= y && wakeLines[i+1]->getY() > y))
            {
                wl1 = wakeLines[i];
                wl2 = wakeLines[i+1];
            }
        }
    }
    double interp = (y-wl1->getY())/(wl2->getY()-wl1->getY());
    double strength = wl1->getStrength()+interp*(wl2->getStrength()-wl1->getStrength());
    return strength;
}

double wake::downwash(Eigen::Matrix<double,1,3> POI, Eigen::MatrixXd pntCloud,Eigen::Matrix<bool,Eigen::Dynamic,1> upperFlag, bool POIflag)
{
//    Eigen::Vector3d V;
    double w;
    Eigen::MatrixXd Xf(pntCloud.rows(),1);
    double phiPOI = 0;
    Eigen::VectorXd dPhi(pntCloud.rows());
    for (int i=0; i<wpanels.size(); i++)
    {
        phiPOI += wpanels[i]->panelPhi(POI);
    }
    if (POIflag)
    {
        phiPOI -= wakeStrength(POI(1));
    }
//    std::cout << "phiPOI = " << phiPOI << std::endl;
    for (int i=0; i<pntCloud.rows(); i++)
    {
        double phiPnt = 0;
        if (upperFlag(i))
        {
            Xf(i) = (pntCloud.row(i)-POI).norm();
        }
        else
        {
            Xf(i) = -(pntCloud.row(i)-POI).norm();
        }
        for (int j=0; j<wpanels.size(); j++)
        {
            phiPnt += wpanels[j]->panelPhi(pntCloud.row(i));
        }
        if (upperFlag(i))
        {
            // Correct for jump in discontinuity;
            phiPnt -= wakeStrength(pntCloud(i,1));
        }
//        meanPhi += phiPnt/pntCloud.rows();
//        std::cout << Xf(i) << "\t" << phiPnt << std::endl;
        dPhi(i) = phiPnt-phiPOI;
    }
//    std::cout << meanPhi << std::endl;
    
    std::ofstream fid;
    fid.open("/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/wakeTLSinputs.txt");
    fid << POI(0) << "\t" << POI(1) << "\t" << POI(2) << "\n";
    fid << pntCloud.rows() << "\n";
    for (int i=0; i<pntCloud.rows(); i++)
    {
        fid << pntCloud(i,0) << "\t" << pntCloud(i,1) << "\t" << pntCloud(i,2) << "\n";
    }
    for (int i=0; i<dPhi.rows(); i++)
    {
        fid << dPhi(i) << "\n";
    }
    fid.close();
    
    Eigen::MatrixXd Xb(0,3),Vb(0,3);
    Eigen::Vector3d V0 = Eigen::Vector3d::Zero();
    Eigen::Matrix<double,1,1> x0;
    x0.setZero();
    chtlsnd weights(x0,Xf,3,Xb,Vb,V0);
//    V(0) = weights.getF().row(0)*dPhi;
//    V(1) = weights.getF().row(1)*dPhi;
//    V(2) = weights.getF().row(2)*dPhi;
    w = weights.getF().row(0)*dPhi;
    return w;
}

Eigen::Vector3d wake::pntVel(Eigen::Matrix<double,1,3> POI, Eigen::MatrixXd pntCloud,Eigen::Matrix<bool,Eigen::Dynamic,1> upperFlag, Eigen::Vector3d Vinf)
{
    //    Eigen::Vector3d V;
    Eigen::Vector3d vel;
    double phiPOI = 0;
    Eigen::VectorXd dPhi(pntCloud.rows());
    for (int i=0; i<wpanels.size(); i++)
    {
        phiPOI += wpanels[i]->panelPhi(POI);
    }
    //    std::cout << "phiPOI = " << phiPOI << std::endl;
    for (int i=0; i<pntCloud.rows(); i++)
    {
        double phiPnt = 0;
        for (int j=0; j<wpanels.size(); j++)
        {
            phiPnt += wpanels[j]->panelPhi(pntCloud.row(i));
        }
        if (POI(2) < 0 && upperFlag(i))
        {
            // Correct for jump in discontinuity;
            phiPnt -= wakeStrength(pntCloud(i,1));
        }
        else if (POI(2) > 0 && !upperFlag(i))
        {
            phiPnt += wakeStrength(pntCloud(i,1));
        }
        //        meanPhi += phiPnt/pntCloud.rows();
        //        std::cout << Xf(i) << "\t" << phiPnt << std::endl;
        dPhi(i) = phiPnt-phiPOI;
        std::cout << pntCloud(i,2) << "\t" << (pntCloud.row(i)-POI).norm() << "\t" << dPhi(i) << std::endl;

    }
    //    std::cout << meanPhi << std::endl;
    
    std::ofstream fid;
    fid.open("/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/wakeTLSinputs.txt");
    fid << POI(0) << "\t" << POI(1) << "\t" << POI(2) << "\n";
    fid << pntCloud.rows() << "\n";
    for (int i=0; i<pntCloud.rows(); i++)
    {
        fid << pntCloud(i,0) << "\t" << pntCloud(i,1) << "\t" << pntCloud(i,2) << "\n";
    }
    for (int i=0; i<dPhi.rows(); i++)
    {
        fid << dPhi(i) << "\n";
    }
    fid.close();
    
    Eigen::MatrixXd Xb(0,3),Vb(0,3);
    Eigen::Vector3d V0 = Eigen::Vector3d::Zero();
    chtlsnd weights(POI,pntCloud,3,Xb,Vb,V0);
    vel(0) = -weights.getF().row(0)*dPhi;
    vel(1) = -weights.getF().row(1)*dPhi;
    vel(2) = -weights.getF().row(2)*dPhi;
    vel += Vinf;
    return vel;
}



void wake::horseshoeTrefftz(double Vinf,double Sref, double &CL, double &CD, Eigen::VectorXd &yLoc, Eigen::VectorXd &Cl, Eigen::VectorXd &Cd)
{
    Eigen::VectorXd Vys = Eigen::VectorXd::Zero(wakeLines.size());
    Eigen::VectorXd Vzs = Eigen::VectorXd::Zero(wakeLines.size());
    yLoc = Eigen::VectorXd::Zero(wakeLines.size());
    Cl = Eigen::VectorXd::Zero(wakeLines.size());
    Cd = Eigen::VectorXd::Zero(wakeLines.size());
    double Vy,Vz,DZT,DYT;
    
    for (int i=0; i<wakeLines.size(); i++)
    {
        for (int j=0; j<wakeLines.size(); j++)
        {
//            if (i != j)
//            {
                wakeLines[j]->horseshoeW(wakeLines[i]->getPMid(), Vy, Vz);
                Vys(i) += Vy;
                Vzs(i) += Vz;
//            }
        }
        yLoc(i) = wakeLines[i]->getY();
        DZT = wakeLines[i]->getP2()(2)-wakeLines[i]->getP1()(2);
        DYT = wakeLines[i]->getP2()(1)-wakeLines[i]->getP1()(1);
        std::cout << yLoc(i) << "\t" << Vys(i) << "\t" << Vzs(i) << std::endl;
        Cd(i) = wakeLines[i]->getStrength()*(DZT*Vys(i)-DYT*Vzs(i))/(Vinf*Vinf*Sref);
        Cl(i) = -2*wakeLines[i]->getStrength()*DYT/(Vinf*Sref);
    }
    CL = Cl.sum();
    CD = Cd.sum();
    
    
}

void wake::sheetTrefftz(double Vinf,double Sref, double &CL, double &CD, Eigen::VectorXd &yLoc, Eigen::VectorXd &Cl, Eigen::VectorXd &Cd)
{
    yLoc = Eigen::VectorXd::Zero(vortexSheets.size());
    Cl = Eigen::VectorXd::Zero(vortexSheets.size());
    Cd = Eigen::VectorXd::Zero(vortexSheets.size());
    double w,gamma;
    
    for (int i=0; i<vortexSheets.size(); i++)
    {
        Eigen::VectorXi verts = vortexSheets[i]->getVerts();
        double dy = (*nodes)(verts(0),1)-(*nodes)(verts(1),1);
        gamma = vortexSheets[i]->getMu();
//        w = sheetDownwash(i);
        w = sheetVradial(i);
        yLoc(i) = vortexSheets[i]->getCenter()(1);
        std::cout << yLoc(i) << "\t" << w << "\t" << gamma << std::endl;
        Cd(i) = -gamma*w*dy/(Vinf*Vinf*Sref);
        Cl(i) = -2*gamma*dy/(Vinf*Sref);
    }
    CL = Cl.sum();
    CD = Cd.sum();
}

double wake::sheetDownwash(int sheetNum)
{
    //    Eigen::Vector3d V;
    Eigen::Vector3d POI = vortexSheets[sheetNum]->getCenter();
    Eigen::Vector3d normal = vortexSheets[sheetNum]->getNormal();
    double w;
    int nPnts = 10;
    double dz = 1;
    double step = 2*dz/(nPnts-1);
    Eigen::MatrixXd Xf(nPnts,1);
    double phiPOI = 0;
    Eigen::VectorXd dPhi(nPnts);
    for (int i=0; i<vortexSheets.size(); i++)
    {
        phiPOI += vortexSheets[i]->panelPhi(POI);
    }
    std::cout << "phiPOI = " << phiPOI << std::endl;
    for (int i=0; i<nPnts; i++)
    {
        double phiPnt = 0;
        double delta = -dz+i*step;
        Eigen::Vector3d pnt = POI+delta*normal;
        Xf(i) = delta;
        for (int j=0; j<vortexSheets.size(); j++)
        {
            phiPnt += vortexSheets[j]->panelPhi(pnt);
        }
        if (delta > 0)
        {
            // Correct for jump in discontinuity;
            phiPnt -= vortexSheets[sheetNum]->getMu();
        }
        std::cout << pnt(0) << "\t" << pnt(1) << "\t" << pnt(2) << "\t" << Xf(i) << "\t" << phiPnt << std::endl;
        dPhi(i) = phiPnt-phiPOI;
    }
    
//    std::ofstream fid;
//    fid.open("/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/wakeTLSinputs.txt");
//    fid << POI(0) << "\t" << POI(1) << "\t" << POI(2) << "\n";
//    fid << pntCloud.rows() << "\n";
//    for (int i=0; i<pntCloud.rows(); i++)
//    {
//        fid << pntCloud(i,0) << "\t" << pntCloud(i,1) << "\t" << pntCloud(i,2) << "\n";
//    }
//    for (int i=0; i<dPhi.rows(); i++)
//    {
//        fid << dPhi(i) << "\n";
//    }
//    fid.close();
    
    Eigen::MatrixXd Xb(0,3),Vb(0,3);
    Eigen::Vector3d V0 = Eigen::Vector3d::Zero();
    Eigen::Matrix<double,1,1> x0;
    x0.setZero();
    chtlsnd weights(x0,Xf,3,Xb,Vb,V0);
    //    V(0) = weights.getF().row(0)*dPhi;
    //    V(1) = weights.getF().row(1)*dPhi;
    //    V(2) = weights.getF().row(2)*dPhi;
    w = weights.getF().row(0)*dPhi;
    return w;
}

double wake::sheetVradial(int sheetNum)
{
    //    Eigen::Vector3d V;
    Eigen::Vector3d pWake = vortexSheets[sheetNum]->getCenter();
    double r;
    double theta = M_PI/12;
    Eigen::Vector3d POI;
    POI(0) = pWake(0);
    if (pWake(1) >= 0)
    {
        r = yMax-pWake(1);
        POI(1) = yMax-r*cos(theta);
        
    }
    else
    {
        r = pWake(1)-yMin;
        POI(1) = yMin+r*cos(theta);
    }
    POI(2) = pWake(2)+r*sin(theta);
    
    double Vr;
    int nPnts = 16;
    if (nPnts % 2 != 0)
    {
        nPnts++; //Make even if odd
    }
    double dz = 0.1;
    double step = 2*dz/(nPnts-1);
    double phiPOI = 0;
    Eigen::VectorXd dPhiy(nPnts);
    Eigen::VectorXd dPhiz(nPnts);
    Eigen::MatrixXd dY(nPnts,1);
    Eigen::MatrixXd dZ(nPnts,1);
    for (int i=0; i<vortexSheets.size(); i++)
    {
        phiPOI += vortexSheets[i]->panelPhi(POI);
    }
    std::cout << "phiPOI = " << phiPOI << std::endl;
    int i=0;
    while (i < nPnts)
    {
        double phiPnt1 = 0;
        double phiPnt2 = 0;
        double delta = -dz+i*step;
        Eigen::Vector3d ydir,zdir;
        ydir << 0,1,0;
        zdir << 0,0,1;
        Eigen::Vector3d pnt1 = POI+delta*ydir;
        Eigen::Vector3d pnt2 = POI+delta*zdir;
        for (int j=0; j<vortexSheets.size(); j++)
        {
            phiPnt1 += vortexSheets[j]->panelPhi(pnt1);
            phiPnt2 += vortexSheets[j]->panelPhi(pnt2);
        }
        if (pnt1(2) > 0)
        {
            // Correct for jump in discontinuity;
            phiPnt1 -= vortexSheets[sheetNum]->getMu();
        }
        if (pnt2(2) > 0)
        {
            phiPnt2 -= vortexSheets[sheetNum]->getMu();
        }
        dPhiy(i) = phiPnt1-phiPOI;
        dPhiz(i) = phiPnt2-phiPOI;
        dY(i) = pnt1(1) - POI(1);
        dZ(i) = pnt2(2) - POI(2);
        std::cout << dY(i) << "\t" << dPhiy(i) << "\t" << dZ(i) << "\t" << dPhiz(i) << std::endl;
        i = i+1;
    }
    
    Eigen::MatrixXd Xb(0,3),Vb(0,3);
    Eigen::Vector3d V0 = Eigen::Vector3d::Zero();
    Eigen::Matrix<double,1,1> x0;
    x0.setZero();
    chtlsnd weightsY(x0,dY,3,Xb,Vb,V0);
    double v = weightsY.getF().row(0)*dPhiy;
    chtlsnd weightsZ(x0,dZ,3,Xb,Vb,V0);
    double w = weightsZ.getF().row(0)*dPhiz;
    Vr = sqrt(pow(v,2)+pow(w,2));
    return Vr;
}

double wake::Vradial(Eigen::Vector3d pWake)
{
    double r;
    double theta = M_PI/8;
//    double deltaZ = 0.5;
    Eigen::Vector3d POI;
    POI(0) = pWake(0);
    if (pWake(1) >= 0)
    {
        r = yMax-pWake(1);
//        if (r < deltaZ)
//        {
//            deltaZ = r;
//        }
//        POI(1) = yMax-sqrt(pow(r,2)-pow(deltaZ,2));
        POI(1) = yMax-r*cos(theta);

    }
    else
    {
        r = pWake(1)-yMin;
//        if (r < deltaZ)
//        {
//            deltaZ = r;
//        }
//        POI(1) = yMin+sqrt(pow(r,2)-pow(deltaZ,2));
        POI(1) = yMin+r*cos(theta);

    }
    POI(2) = pWake(2)+r*sin(theta);
//    POI(2) = pWake(2)+deltaZ;
    
    double Vr;
    int nPnts = 16;
    if (nPnts % 2 != 0)
    {
        nPnts++; //Make even if odd
    }
    double dz = 0.05;
    double step = 2*dz/(nPnts-1);
    double phiPOI = 0;
    Eigen::VectorXd dPhiy(nPnts);
    Eigen::VectorXd dPhiz(nPnts);
    Eigen::MatrixXd dY(nPnts,1);
    Eigen::MatrixXd dZ(nPnts,1);
    for (int i=0; i<wpanels.size(); i++)
    {
        phiPOI += wpanels[i]->panelPhi(POI);
    }
//    std::cout << "phiPOI = " << phiPOI << std::endl;
    int i=0;
    while (i < nPnts)
    {
        double phiPnt1 = 0;
        double phiPnt2 = 0;
        double delta = -dz+i*step;
        Eigen::Vector3d ydir,zdir;
        ydir << 0,1,0;
        zdir << 0,0,1;
        Eigen::Vector3d pnt1 = POI+delta*ydir;
        Eigen::Vector3d pnt2 = POI+delta*zdir;
        for (int j=0; j<wpanels.size(); j++)
        {
            phiPnt1 += wpanels[j]->panelPhi(pnt1);
            phiPnt2 += wpanels[j]->panelPhi(pnt2);
        }
        if (pnt1(2) > 0)
        {
            // Correct for jump in discontinuity;
            phiPnt1 -= wakeStrength(pWake(1));
        }
        if (pnt2(2) > 0)
        {
            phiPnt2 -= wakeStrength(pWake(1));
        }
        dPhiy(i) = phiPnt1-phiPOI;
        dPhiz(i) = phiPnt2-phiPOI;
        dY(i) = pnt1(1) - POI(1);
        dZ(i) = pnt2(2) - POI(2);
        if (r < .31)
        {
            std::cout << pnt1(0) << "\t" << pnt1(1) << "\t" << pnt1(2) << "\t" << phiPnt1 << "\t" << pnt2(0) << "\t" << pnt2(1) << "\t" << pnt2(2) << "\t" << phiPnt2 << std::endl;
        }
//        std::cout << dY(i) << "\t" << dPhiy(i) << "\t" << dZ(i) << "\t" << dPhiz(i) << std::endl;
        i = i+1;
    }
    
    Eigen::MatrixXd Xb(0,3),Vb(0,3);
    Eigen::Vector3d V0 = Eigen::Vector3d::Zero();
    Eigen::Matrix<double,1,1> x0;
    x0.setZero();
    chtlsnd weightsY(x0,dY,3,Xb,Vb,V0);
    double v = weightsY.getF().row(0)*dPhiy;
    chtlsnd weightsZ(x0,dZ,3,Xb,Vb,V0);
    double w = weightsZ.getF().row(0)*dPhiz;
    Vr = sqrt(pow(v,2)+pow(w,2));
    return Vr;
}

Eigen::Vector3d wake::pntInWake(double x, double y)
{
    Eigen::Vector3d p1,p2,tvec,pnt,out,pntInWake;
    double t,scale;
    pntInWake.setZero();
    for (int i=0; i<vortexSheets.size(); i++)
    {
        Eigen::VectorXi verts = vortexSheets[i]->getVerts();
        p1 = nodes->row(verts(0));
        p2 = nodes->row(verts(1));
        if (p1(1) >= y && p2(1) < y)
        {
            t = (y-p1(1))/(p2(1)-p1(1));
            tvec = (p2-p1);
            pnt = p1+t*tvec;
            out = vortexSheets[i]->getNormal().cross(tvec);
            scale = (x-pnt(0))/out(0);
            pntInWake = pnt+scale*out;
            break;
        }
    }
    return pntInWake;
    
}