//
//  bodyPanel.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "bodyPanel.h"

void bodyPanel::setNeighbors(panelOctree *oct, short normalMax)
{
    short absMax = normalMax;
    if (upper || lower)
    {
        absMax = verts.size()+1;
    }
    else
    {
        absMax = normalMax;
    }
    node<panel>* currentNode = oct->findNodeContainingMember(this);
    node<panel>* exception = NULL;
    
    while (neighbors.size() < absMax)
    {
        scanForNeighbors(currentNode,exception);
        if (currentNode == oct->getRootNode())
        {
            break;
        }
        else
        {
            exception = currentNode;
            currentNode = currentNode->getParent();
        }
    }
}

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

std::vector<bodyPanel*> bodyPanel::getBodyNeighbors()
{
    std::vector<bodyPanel*> bodyNeighbors;
    for (int i=0; i<neighbors.size(); i++)
    {
        if (neighbors[i]->getID() < 10000)
        {
            bodyPanel* p = dynamic_cast<bodyPanel*>(neighbors[i]);
            bodyNeighbors.push_back(p);
        }
    }
    return bodyNeighbors;
}

std::vector<bodyPanel*> bodyPanel::gatherNeighbors(int nPanels)
{
    std::vector<bodyPanel*> cluster;
    unsigned long oldSize = cluster.size();
    cluster.push_back(this);
    bool tipFlag = wingTipTest(this);
    while (cluster.size() < nPanels+1)
    {
        std::vector<bodyPanel*> toAdd;
        for (unsigned long i=oldSize; i<cluster.size(); i++)
        {
            std::vector<bodyPanel*> temp = cluster[i]->getBodyNeighbors();
            for (int j=0; j<temp.size(); j++)
            {
                if (tipFlag)
                {
                    if (clusterTest(temp[j], 0.1, cluster))
                    {
                        // If panel is on wing tip, only include panels also on wing tip
                        toAdd.push_back(temp[j]);
                    }
                }
                else if (nearTrailingEdge())
                {
                    if (clusterTest(temp[j], M_PI/2.1, cluster))
                    {
                        toAdd.push_back(temp[j]);
                    }
                }
                else if (clusterTest(temp[j], 5*M_PI/6, cluster))
                {
                    // Do not include panels on other side of discontinuity (i.e. wake), panels already in cluster, this panel, or wake panels
                    toAdd.push_back(temp[j]);
                }
            }
        }
        assert(cluster.size() != oldSize);
        oldSize = cluster.size();
        for (int i=0; i<toAdd.size(); i++)
        {
            cluster.push_back(toAdd[i]);
            if (cluster.size() == nPanels+1)
            {
                break;
            }
        }
    }
    cluster.erase(cluster.begin());
    bool flag = false;
    if (wingTipTest(this) && flag == false)
    {
        flag = true;
        std::ofstream fid;
        fid.open("/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/ClusterCheck_lower.txt");
        fid << center(0) << "\t" << center(1) << "\t" << center(2) << "\n";
        for (int i=0; i<cluster.size(); i++)
        {
            fid << cluster[i]->getCenter()(0) << "\t" << cluster[i]->getCenter()(1) << "\t" << cluster[i]->getCenter()(2) << "\n";
        }
        fid.close();
    }
    return cluster;
}

void bodyPanel::computeVelocity()
{
    int TSorder = 3;
    int obs,dim;
    std::vector<bodyPanel*> cluster;
    Eigen::MatrixXd Xf,Xb,Vb;
    Eigen::VectorXd df;
    if (wingTipTest(this))
    {
        Eigen::Vector2d X0 = Eigen::Vector2d::Zero();
        Eigen::Vector2d V0 = Eigen::Vector2d::Zero();
        dim = 2;
        obs = chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder)); // Binomial Coefficient
        obs += 5;
        cluster = gatherNeighbors(obs);
        Xf.resize(obs,3);
        Xb = Eigen::MatrixXd::Zero(0,dim);
        Vb = Eigen::MatrixXd::Zero(0,dim);
        df.resize(obs);
        for (int i=0; i<obs; i++)
        {
            Xf.row(i) = global2local(cluster[i]->getCenter(),true);
            df(i) = cluster[i]->getPotential()-potential;
        }
        Eigen::MatrixXd xLocal = Xf.block(0,0,obs,2);
        chtlsnd tipV(X0,xLocal,TSorder,Xb,Vb,V0);
        Eigen::Vector3d vLocal;
        vLocal(0) = tipV.getF().row(0)*df;
        vLocal(1) = tipV.getF().row(1)*df;
        vLocal(2) = 0;
        velocity = local2global(vLocal,false);
    }
    else
    {
        dim = 3;
        obs = chtlsnd::factorial(TSorder+dim)/(chtlsnd::factorial(dim)*chtlsnd::factorial(TSorder)); // Binomial Coefficient
        obs += 5;
        cluster = gatherNeighbors(obs);
        Xf.resize(obs,dim);
        Xb.resize(obs,dim);
        Vb.resize(obs,dim);
        df.resize(obs);
        for (int i=0; i<obs; i++)
        {
            Xf.row(i) = cluster[i]->getCenter();
            Vb.row(i) = cluster[i]->getBezNormal();
            df(i) = cluster[i]->getPotential()-potential;
        }
        Xb = Xf;
        chtlsnd vWeights(center,Xf,TSorder,Xb,Vb,bezNormal);
        velocity(0) = vWeights.getF().row(0)*df;
        velocity(1) = vWeights.getF().row(1)*df;
        velocity(2) = vWeights.getF().row(2)*df;
    }
}

void bodyPanel::computeCp(double Vinf)
{
    Cp = 1-pow(velocity.norm()/Vinf,2);
}

bool bodyPanel::clusterTest(bodyPanel* other,double angle, const std::vector<bodyPanel*> &cluster)
{
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

bool bodyPanel::wingTipTest(bodyPanel* p)
{
    if (!p->isLiftSurf())
    {
        return false;
    }
    else
    {
        // Catch panels with a wingtip normal component in y direction corresponding to max sweep and dihedral of 45 and 15 degrees, respectively.
        return std::abs(p->getNormal()(1)) >= .7;
    }
    
}

bool bodyPanel::nearTrailingEdge()
{
    // Returns true if panel is on trailing edge or a neighbor is on trailing edge
    if (lsFlag && (upper || lower))
    {
        return true;
    }
    else if (lsFlag)
    {
        std::vector<bodyPanel*> neighbs = getBodyNeighbors();
        for (int i=0; i<neighbs.size(); i++)
        {
            if (neighbs[i]->isUpper() || neighbs[i]->isLower())
            {
                return true;
            }
        }
    }
    return false;
}

void bodyPanel::tipVelocity()
{
    // Computes velocity for panels on flat wing tip by doing a 2D Taylor Series Least Squares in local coordinate frame.
    
}
