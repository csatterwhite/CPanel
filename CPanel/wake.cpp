//
//  wake.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wake.h"
#include "wakePanel.h"
#include "wakeLine.h"
#include "edge.h"
#include "cpNode.h"


wake::~wake()
{
    for (int i=0; i<wakeLines.size(); i++)
    {
        delete wakeLines[i];
    }
    for (int i=0; i<vortexSheets.size(); i++)
    {
        delete vortexSheets[i];
    }
}

void wake::addPanel(wakePanel* wPan)
{
    Eigen::Vector3d pnt;
    std::vector<cpNode*> nodes = wPan->getNodes();
    if (wpanels.size() == 0)
    {
        // Initialize wake dimensions with first panel
        pnt = nodes[0]->getPnt();
        yMax = pnt(1);
        yMin = yMax;
        x0 = pnt(0);
        xf = x0;
        z0 = pnt(2);
        zf = z0;
        normal = wPan->getNormal();
    }
    
    for (int i=0; i<nodes.size(); i++)
    {
        pnt = nodes[i]->getPnt();
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
        }
        else if (pnt(0) < x0)
        {
            x0 = pnt(0);
        }
        if (pnt(2) > zf)
        {
            zf = pnt(2);
        }
        else if (pnt(2) < z0)
        {
            z0 = pnt(2);
        }
    }
    wpanels.push_back(wPan);
}

bool wake::isSameWake(wake* other)
{
    if (other == this)
    {
        return false;
    }
    
//    Eigen::Vector3d p1 = wpanels[0]->getCenter();
//    Eigen::Vector3d p2 = other->getPanels()[0]->getCenter();
//    
//    Eigen::Vector3d vec = p2-p1;
//    
//    double dot = vec.dot(normal)/(vec.norm());
//    if (std::abs(dot) < pow(10,-10) && (yMin == other->getYMax() || yMax == other->getYMin()))
//    {
//        // Normal vector and vector connecting two points are orthogonal
//        return true;
//    }
    
    double eps = pow(10, -4);
//    if (other->getX0() == x0 && other->getZ0() == z0 && other->getXf() == xf && other->getZf() == zf)
    if (std::abs(other->getX0() - x0) < eps && std::abs(other->getZ0() - z0) < eps && std::abs(other->getXf() - xf) < eps && std::abs(other->getZf() - zf) < eps)
    {
        return true;
    }
    
    return false;
}

void wake::mergeWake(wake *other)
{
    std::vector<wakePanel*> pans = other->getPanels();
    wakePanel* w;
    for (int i=0; i<pans.size(); i++)
    {
        w = pans[i];
        wpanels.push_back(w);
        w->setParentWake(this);
    }
    
    std::vector<wakeLine*> otherLines = other->getWakeLines();
    wakeLine* wLine;
    for (int i=0; i<otherLines.size(); i++)
    {
        wLine = new wakeLine(*otherLines[i]);
        addWakeLine(wLine);
    }
    
//    ///////
//    std::cout << std::endl;
//    for (int i=0; i<wakeLines.size(); i++)
//    {
//        std::cout << wakeLines[i]->getPMid()(0) << "," << wakeLines[i]->getPMid()(1) << "," << wakeLines[i]->getPMid()(2) << ";" << std::endl;
//    }
//    ////
    
    if (other->getYMin() < yMin)
    {
        yMin = other->getYMin();
    }
    if (other->getYMax() > yMax)
    {
        yMax = other->getYMax();
    }
}

void wake::addTEPanel(wakePanel* p)
{
    TEpanels.push_back(p);
}

void wake::addWakeLine(wakeLine* wl)
{
    wakeLines.push_back(wl);
    std::sort(wakeLines.begin(),wakeLines.end(),[](wakeLine* w1, wakeLine* w2) {return w1->getY() < w2->getY();});
}

void wake::trefftzPlane(double Vinf,double Sref)
{
    int nPnts = 60;
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

    double xTrefftz = x0+2*(xf-x0)/3;
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

Eigen::Vector3d wake::lambVectorInt(const Eigen::Vector3d &Vinf,Eigen::VectorXd &yLoc)
{
    // Sort by y position
    std::sort(TEpanels.begin(), TEpanels.end(), [](wakePanel* w1, wakePanel* w2) {return w1->getCenter()(1) < w2->getCenter()(1);});
    yLoc.resize(TEpanels.size()+2);
    edge* TE = TEpanels[0]->getTE();
    
    if (TE->getN1()->getPnt()(1) > TE->getN2()->getPnt()(1))
    {
        TE->flipDir();
    }
    
    int i = 1;
    Eigen::Vector3d vel,circ;
    Eigen::MatrixXd sectForces = Eigen::MatrixXd::Zero(TEpanels.size()+2,3);
    while (TE != nullptr)
    {
        yLoc(i) = TE->getMidPoint()(1);
        vel = TE->edgeVelocity(Vinf);
        circ = TE->TEgamma();
        sectForces.row(i) = vel.cross(circ);
        TE = TE->nextTE();
        i++;
    }
    
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    Eigen::Vector3d sectF1,sectF2;
    i = 0;
    
    while (i < sectForces.rows()-1)
    {
        sectF1 = sectForces.row(i);
        sectF2 = sectForces.row(i+1);
        F = F + 0.5*(yLoc(i+1)-yLoc(i))*(sectF1+sectF2);
        i++;
    }
    
    return F;
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
            phiPnt1 += wakeStrength(pWake(1));
        }
        if (pnt2(2) > 0)
        {
            phiPnt2 += wakeStrength(pWake(1));
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
    std::vector<edge*> edges;
    pntInWake.setZero();
    for (int i=0; i<wpanels.size(); i++)
    {
        if (wpanels[i]->isTEpanel())
        {
            std::vector<edge*> edges = wpanels[i]->getUpper()->getEdges();
            for (int j=0; j<edges.size(); j++)
            {
                if (edges[j]->isTE())
                {
                    p1 = edges[j]->getNodes()[0]->getPnt();
                    p2 = edges[j]->getNodes()[1]->getPnt();
                    if ((p1(1) <= y && p2(1) >= y) || (p1(1) >= y && p2(1) <= y))
                    {
                        t = (y-p1(1))/(p2(1)-p1(1));
                        tvec = (p2-p1);
                        pnt = p1+t*tvec;
                        out = wpanels[i]->getNormal().cross(tvec);
                        if (out(0) < 0)
                        {
                            out *= -1; // Flip sign if p1 and p2 were out of order
                        }
                        scale = (x-pnt(0))/out(0);
                        pntInWake = pnt+scale*out;
                        return pntInWake;
                    }
                }
            }
        }
    }
    return pntInWake;
}