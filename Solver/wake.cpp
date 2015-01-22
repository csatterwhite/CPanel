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
    int nPnts = 30;
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