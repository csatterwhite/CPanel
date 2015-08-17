//
//  bodyPanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "bodyPanel.h"
#include "edge.h"
#include "cpNode.h"
#include "surface.h"

bodyPanel::bodyPanel(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm,surface* parentSurf, int surfID) : panel(nodes,pEdges,bezNorm,surfID), upper(false), lower(false), streamFlag(false), parentSurf(parentSurf), TEpanel(false), TSorder(3), tipFlag(false)
{
    for (int i=0; i<pEdges.size(); i++)
    {
        pEdges[i]->addBodyPan(this);
    }
    for (int i=0; i<nodes.size(); i++)
    {
        nodes[i]->addBodyPanel(this);
    }
}

//bodyPanel::bodyPanel(const bodyPanel &copy) : panel(copy), sourceStrength(copy.sourceStrength), tipFlag(copy.tipFlag) {}

void bodyPanel::addNeighbor(bodyPanel* p)
{
    neighbors.push_back(p);
}

void bodyPanel::setUpper() {upper = true;}
void bodyPanel::setLower() {lower = true;}

void bodyPanel::setTEpanel(edge* trailingEdge)
{
    TEpanel = true;
    parentSurf->setLSflag();
    TE = trailingEdge;
}

void bodyPanel::setIndex(int i) {index = i;}

void bodyPanel::setTipFlag()
{
    if (!parentSurf->isLiftingSurf())
    {
        tipFlag = false;
        return;
    }
    
    int count = 0;
    double angle;
    Eigen::Vector3d nNormal;
    for (int i=0; i<neighbors.size(); i++)
    {
        nNormal = neighbors[i]->getNormal();
        double dot = normal.dot(nNormal)/(normal.norm()*nNormal.norm());
        if (dot > 1)
        {
            dot = 1;
        }
        else if (dot < -1)
        {
            dot = -1;
        }
        
        angle = acos(dot);
        
        if (angle < pow(10,-8) && asin(nNormal(0)) > -M_PI/12)
        {
            // asin(nNormal(0)) is the angle normal vector makes with xy plane.  Catches bug where some leading edge panels were being flagged as on the tip patch.
            count++;
        }
    }
    if (count >= 2)
    {
        tipFlag = true;
    }
    
}

void bodyPanel::setSigma(Eigen::Vector3d Vinf, double Vnorm)
{
    sourceStrength = (-Vinf.dot(normal)+Vnorm);
}

void bodyPanel::setMu(double dubStrength)
{
    doubletStrength = dubStrength;
}

void bodyPanel::setStreamFlag()
{
    streamFlag = true;
}

double bodyPanel::panelPhi(const Eigen::Vector3d &POI)
{
    double phiSrc,phiDub;
    phiSrc = 0;
    phiDub = 0;
    
    panelPhiInf(POI,phiSrc,phiDub);
    return -sourceStrength*phiSrc-doubletStrength*phiDub;
}

Eigen::Vector3d bodyPanel::panelV(const Eigen::Vector3d &POI)
{
    Eigen::Vector3d vSrc,vDub;
    
    panelVInf(POI,vSrc,vDub);
    
    return sourceStrength*vSrc+doubletStrength*vDub;
}

void bodyPanel::panelPhiInf(const Eigen::Vector3d &POI, double &phiSrc,double &phiDub)
{
    Eigen::Vector3d pjk = POI-center;
    bool itselfFlag = false;
    if (pjk.norm() < pow(10,-10))
    {
        itselfFlag = true;
    }
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
        for (int i=0; i<nodes.size(); i++)
        {
            Eigen::Vector3d p1;
            Eigen::Vector3d p2;
            if (i!=nodes.size()-1)
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[i+1]->getPnt();
            }
            else
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[0]->getPnt();
            }
            a = POI-p1;
            b = POI-p2;
            s = p2-p1;
            Al = local.row(2).dot(s.cross(a));
            if (!itselfFlag)
            {
                phiV = vortexPhi(PN,Al,a,b,s,local.row(0),local.row(1),local.row(2));
                phiDub += phiV;
            }
            phiSrc += srcSidePhi(PN,Al,phiV,a,b,s);
            
        }
        phiSrc /= (4*M_PI);
        if (!itselfFlag)
        {
            phiDub /= (4*M_PI);
        }
        else
        {
            phiDub = -0.5;
        }
    }
}

void bodyPanel::panelVInf(const Eigen::Vector3d &POI, Eigen::Vector3d &vSrc,Eigen::Vector3d &vDub)
{
    
    // VSAero source and doublet velocity influence formulation
    Eigen::Vector3d pjk = POI-center;
    Eigen::Matrix3d local = getLocalSys();
    double PN = pjk.dot(local.row(2));
    if (pjk.norm()/longSide > 5)
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
        for (int i=0; i<nodes.size(); i++)
        {
            if (i!=nodes.size()-1)
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[i+1]->getPnt();
            }
            else
            {
                p1 = nodes[i]->getPnt();
                p2 = nodes[0]->getPnt();
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

void bodyPanel::setCluster()
{
    int dim;
    int buffer = 10;
    if (tipFlag)
    {
        dim = 2;
    }
    else
    {
        dim = 3;
    }
    
    double nObs = chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder))-1+buffer; // Binomial Coefficient
    double nPanels = std::ceil(nObs/2);
    int oldSize = (int)cluster.size();
    cluster.push_back(this);
    bool upFlag = upper;
    bool lowFlag = lower;
    while (cluster.size() < nPanels)
    {
        std::vector<bodyPanel*> toAdd;
        for (unsigned long i=oldSize; i<cluster.size(); i++)
        {
            std::vector<bodyPanel*> temp = cluster[i]->getNeighbors();
            for (int j=0; j<temp.size(); j++)
            {
                if (clusterTest(temp[j], 5*M_PI/6,upFlag,lowFlag))
                {
                    // Do not include panels on other side of discontinuity (i.e. wake), panels already in cluster, or this panel
                    if (temp[j]->isUpper())
                    {
                        upFlag = true;
                    }
                    else if (temp[j]->isLower())
                    {
                        lowFlag = true;
                    }
                    if (std::find(toAdd.begin(), toAdd.end(), temp[j]) == toAdd.end())
                    {
                        toAdd.push_back(temp[j]);
                    }
                }
            }
        }
        
        if (toAdd.size() == 0 && cluster.size() < nPanels)
        {
            // No valid panels could be found so decrease order of taylor series
            while (cluster.size() < nPanels+1 && TSorder >= 1)
            {
                if (buffer > 0)
                {
                    buffer -= 1;
                }
                else
                {
                    TSorder -= 1;
                }
                
                // Do Something to fix error when cluster can't be found.
                
                nObs = chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder)) + buffer; // Binomial Coefficient
                nPanels = ceil(nObs/2);
            }
            
            cluster.erase(cluster.begin()+nPanels+1,cluster.end());
            break;
        }
        else
        {
            oldSize = (int)cluster.size();
            for (int i=0; i<toAdd.size(); i++)
            {
                cluster.push_back(toAdd[i]);
                if (cluster.size() == nPanels+1)
                {
                    break;
                }
            }
        }
    }
}


Eigen::Vector3d bodyPanel::pntVelocity(const Eigen::Vector3d &pnt,double pntPotential, double PG, const Eigen::Vector3d &Vinf)
{
    Eigen::Vector3d vel;
    
    if (cluster.size() == 0)
    {
        setCluster();
    }
    
    if (tipFlag)
    {
        vel = velocity2D(pnt,pntPotential);
    }
    else
    {
        vel = velocity3D(pnt,pntPotential);
        if (vel != vel)
        {
            // NaN can be returned if supporting data was all on a flat patch (essentially 2D)
            vel = velocity2D(pnt,pntPotential);
        }
    }
    
    if (vel != vel)
    {
        // Last resort approximation if CHTLS not able to compute velocity. From Kinney (CBAERO).
        vel = (bezNormal.cross(Vinf)).cross(bezNormal);
    }
    
    vel(0) /= PG; // Prandlt Glauert Correction
    return vel;
}

Eigen::Vector3d bodyPanel::velocity2D(const Eigen::Vector3d &pnt,double pntPotential)
{
    int dim = 2;
    std::vector<bodyPanel*> clust = cluster;
    if (pnt == center)
    {
        clust.erase(clust.begin());
    }
    Eigen::Vector3d vel;
    Eigen::MatrixXd Xf,Xb,Vb;
    Eigen::VectorXd df;
    
    Eigen::Vector3d pntLocal = global2local(pnt,true);
    Eigen::Vector2d X0 = pntLocal.head(2);
    Eigen::Vector2d V0 = Eigen::Vector2d::Zero();
    Xf.resize(clust.size(),3);
    Xb = Eigen::MatrixXd::Zero(0,dim);
    Vb = Eigen::MatrixXd::Zero(0,dim);
    df.resize(clust.size());
    for (int i=0; i<clust.size(); i++)
    {
        Xf.row(i) = global2local(clust[i]->getCenter(),true);
        df(i) = clust[i]->getPotential()-pntPotential;
    }
    Eigen::MatrixXd xLocal = Xf.block(0,0,clust.size(),2);
    chtlsnd tipV(X0,xLocal,TSorder,Xb,Vb,V0);
    Eigen::Vector3d vLocal;
    vLocal(0) = tipV.getF().row(0)*df;
    vLocal(1) = tipV.getF().row(1)*df;
    vLocal(2) = 0;
    vel = local2global(vLocal,false);
    
    return vel;
}

Eigen::Vector3d bodyPanel::velocity3D(const Eigen::Vector3d &pnt,double pntPotential)
{
    int dim = 3;
    std::vector<bodyPanel*> clust = cluster;
    if (pnt == center)
    {
        clust.erase(clust.begin());
    }
    Eigen::Vector3d vel;
    Eigen::MatrixXd Xf,Xb,Vb;
    Eigen::VectorXd df;
    
    Xf.resize(clust.size(),dim);
    Xb.resize(clust.size(),dim);
    Vb.resize(clust.size(),dim);
    df.resize(clust.size());
    for (int i=0; i<clust.size(); i++)
    {
        Xf.row(i) = clust[i]->getCenter();
        Vb.row(i) = clust[i]->getBezNormal();
        df(i) = clust[i]->getPotential()-pntPotential;
    }
    Xb = Xf;
    chtlsnd vWeights(pnt,Xf,TSorder,Xb,Vb,Eigen::Vector3d::Zero());
    
    vel(0) = vWeights.getF().row(0)*df;
    vel(1) = vWeights.getF().row(1)*df;
    vel(2) = vWeights.getF().row(2)*df;
    
    return vel;
}

void bodyPanel::computeCp(double Vinf)
{
    Cp = (1-pow(velocity.norm()/Vinf,2));
}

void bodyPanel::computeVelocity(double PG, const Eigen::Vector3d &Vinf)
{
    velocity = pntVelocity(center,potential,PG,Vinf);
}

Eigen::Vector3d bodyPanel::computeMoments(const Eigen::Vector3d &cg)
{
    Eigen::Vector3d r = center-cg;
    Eigen::Vector3d F = -Cp*bezNormal*area;
    return r.cross(F);
}

bool bodyPanel::clusterTest(bodyPanel* other,double angle,bool upFlag,bool lowFlag)
{
    if ((upFlag && other->isLower()) || (lowFlag && other->isUpper()))
    {
        return false;
    }
    
    if (tipFlag != other->isTipPan())
    {
        return false;
    }
    
    double dot = other->getNormal().dot(normal);
    // Floating point error can cause panels with the same normal vector to result in a dot product greater than one, causing acos to return nan.
    if (dot > 1)
    {
        dot = 1;
    }
    else if (dot < -1)
    {
        dot = -1;
    }
    return (acos(dot) < angle && std::find(cluster.begin(),cluster.end(),other)==cluster.end() && other != this);
}

void bodyPanel::setLSflag()
{
    parentSurf->setLSflag();
}

bool bodyPanel::nearTrailingEdge()
{
    // Returns true if panel is on trailing edge or a neighbor is on trailing edge
    if (parentSurf->isLiftingSurf())
    {
        if (upper || lower)
        {
            return true;
        }
        else
        {
            std::vector<bodyPanel*> neighbs = getNeighbors();
            for (int i=0; i<neighbs.size(); i++)
            {
                if (neighbs[i]->isUpper() || neighbs[i]->isLower())
                {
                    return true;
                }
            }
        }
    }
    return false;
}

double bodyPanel::dist2Pan(bodyPanel* other)
{
    return (other->getCenter()-center).norm();
}

std::vector<bodyPanel*> bodyPanel::getRelatedPanels()
{
    std::vector<bodyPanel*> pans;
    pans.push_back(this);
    std::vector<bodyPanel*> nodePans;
    
    for (int i=0; i<nodes.size(); i++)
    {
        nodePans = nodes[i]->getBodyPans();
        for (int j=0; j<nodePans.size(); j++)
        {
            if (std::find(pans.begin(),pans.end(),nodePans[j]) == pans.end())
            {
                pans.push_back(nodePans[j]);
            }
            
        }
    }
    return pans;
}


