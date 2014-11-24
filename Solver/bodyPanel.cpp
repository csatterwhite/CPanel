//
//  bodyPanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "bodyPanel.h"

double bodyPanel::panelPhi(const Eigen::Vector3d &POI)
{
    double mu,sigma,phiSrc,phiDub;
    
    mu = doubletStrength;
    sigma = sourceStrength;
    
    panelPhiInf(POI,phiSrc,phiDub);
    return sigma*phiSrc+mu*phiDub;
    
}

Eigen::Vector3d bodyPanel::panelV(const Eigen::Vector3d &POI)
{
    double mu,sigma;
    Eigen::Vector3d vSrc,vDub;
    
    mu = doubletStrength;
    sigma = sourceStrength;
    panelVInf(POI,vSrc,vDub);
    
    return sigma*vSrc+mu*vDub;
}

void bodyPanel::panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub)
{
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm()/longSide > 5)
    {
        phiSrc = pntSrcPhi(pjk.norm());
        phiDub = pntDubPhi(PN,pjk.norm());
    }
    else
    {
        double Al,phiV;
        Eigen::Vector3d a,b,s;
        for (int i=0; i<verts.size(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            if (i!=verts.size()-1)
            {
                p1 = nodes->row(verts(i));
                p2 = nodes->row(verts(i+1));
            }
            else
            {
                p1 = nodes->row(verts(i));
                p2 = nodes->row(verts(0));
            }
            a = POI-p1;
            b = POI-p2;
            s = p2-p1;
            Al = local.row(2).dot(s.cross(a));
            phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
            phiSrc += srcSidePhi(PN,Al,phiV,a,b,s);
            phiDub += phiV;
        }
        phiSrc /= (4*M_PI);
        phiDub /= (4*M_PI);
    }
}

void bodyPanel::panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub)
{
    
    // VSAero source and doublet velocity influence formulation
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm() < .00001)
    {
        vSrc = 0.5*normal;
        vDub << 0,0,0;
        return;
    }
    else if (pjk.norm()/longSide > 5)
    {
        vSrc = pntSrcV(pjk);
        vDub = pntDubV(local.row(2),pjk);
    }
    else
    {
        Eigen::Vector3d p1,p2,a,b,s,l,m,n;
        l = local.row(0);
        m = local.row(1);
        n = local.row(2);
        pjk = POI-center;
        double Al;
        for (int i=0; i<verts.size(); i++)
        {
            if (i!=verts.size()-1)
            {
                p1 = nodes->row(verts(i));
                p2 = nodes->row(verts(i+1));
            }
            else
            {
                p1 = nodes->row(verts(i));
                p2 = nodes->row(verts(0));
            }
            a = POI-p1;
            b = POI-p2;
            s = a-b;
            Al = n.dot(s.cross(a));
            
            vDub += vortexV(a,b,s);
            vSrc += srcSideV(PN,Al,a,b,s,l,m,n);
        }
        vDub /= (4*M_PI);
        vSrc /= (4*M_PI);
    }
}

double bodyPanel::srcSidePhi(const double &PN,const double &Al, const double &phiV,const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s)
{
    double A,B,S;
    A = a.norm();
    B = b.norm();
    S = s.norm();
    double GL = 0;
    if (std::abs(A+B-S) > 0 && S > 0)
    {
    	GL = 1/S*log(std::abs((A+B+S)/(A+B-S)));
    }
    return (Al*GL-PN*phiV);
}

Eigen::Vector3d bodyPanel::srcSideV(const double &PN,const double &Al, const Eigen::Vector3d &a,const Eigen::Vector3d &b, const Eigen::Vector3d &s,const Eigen::Vector3d &l,const Eigen::Vector3d &m,const Eigen::Vector3d &n)
{
    double A,B,S;
    A = a.norm();
    B = b.norm();
    S = s.norm();
    double GL = 0;
    if (std::abs(A+B-S) > 0 && S > 0)
    {
        GL = 1/S*log(std::abs((A+B+S)/(A+B-S)));
    }
    double CJK = vortexPhi(PN,Al,a,b,s,l,m,n);
    return (GL*(s.dot(m)*l-s.dot(l)*m)+CJK*n);
}

inline double bodyPanel::pntSrcPhi(const double &PJK)
{
    return area/(4*M_PI*PJK);
}

inline Eigen::Vector3d bodyPanel::pntSrcV(const Eigen::Vector3d &pjk)
{
    return area*pjk/(4*M_PI*pow(pjk.norm(),3));
}

