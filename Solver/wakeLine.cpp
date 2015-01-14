//
//  wakeLine.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakeLine.h"

void wakeLine::setDimensions()
{
    std::vector<Eigen::Vector3d> TEverts;
    Eigen::MatrixXd* nodes = upper->getNodes();
    Eigen::VectorXi uVerts = upper->getVerts();
    Eigen::VectorXi lVerts = lower->getVerts();
    for (int i=0; i<upper->getVerts().size(); i++)
    {
        for (int j=0; j<lower->getVerts().size(); j++)
        {
            if (nodes->row(uVerts(i)) == nodes->row(lVerts(j)))
            {
                TEverts.push_back(nodes->row(uVerts(i)));
            }
        }
    }
    std::sort(TEverts.begin(), TEverts.end(), [](Eigen::Vector3d p1, Eigen::Vector3d p2) {return p1(1)<p2(1);});
    p1 = TEverts[0];
    p2 = TEverts[1];
    pMid = (p1+p2)/2;
}

void wakeLine::surveyPnts(double y,double distDownStream, double x0Wake, double z0Wake, Eigen::Vector3d &trefftzPnt, Eigen::MatrixXd &survPnts, Eigen::Matrix<bool,Eigen::Dynamic,1> &upperFlag)
{
//    double R = (p2-p1).norm(); // Radius of sphere containing survey points
    
    if (survPnts.rows() % 2 != 0)
    {
        survPnts.resize(survPnts.rows()+1,survPnts.cols());
        upperFlag.resize(survPnts.rows()+1, 1);
    }
    double R = 0.5;
    double step = 2*R/survPnts.rows();
    
    
    Eigen::Vector3d pTE; // TE point at y location
    Eigen::Vector3d v = (p2-p1);
    double t = (y-p1(1))/(p2(1)-p1(1));
    pTE = p1+v*t;
    
    double dist = distDownStream - sqrt(pow(pTE(0)-x0Wake,2)+pow(pTE(2)-z0Wake,2));
    
    Eigen::Vector3d normT; // Normal vector of Trefftz plane
    normT(0) = normal(2);
    normT(1) = normal(1);
    normT(2) = -normal(0);
    trefftzPnt = pTE+dist*normT; // Point in Trefftz plane on wake sheet
    Eigen::Vector3d pNew,pVec;
    double r,theta,phi;
    int i = 0;
    while (i < survPnts.rows())
    {
        survPnts.row(i) = trefftzPnt+step*(i+1)*normal;
        survPnts.row(i+1) = trefftzPnt-step*(i+1)*normal;
        upperFlag(i) = true;
        upperFlag(i+1) = false;
        i = i+2;

//        r = R*((double) rand()/(RAND_MAX))+0.1*R;
//        theta = 2*M_PI*((double) rand()/(RAND_MAX));
//        phi = M_PI*((double) rand()/(RAND_MAX));
//        while (phi > M_PI/3 && phi < 2*M_PI/3)
//        {
//            phi = M_PI*((double) rand()/(RAND_MAX));
//        }
//        pNew(0) = r*cos(theta)*sin(phi);
//        pNew(1) = r*sin(theta)*sin(phi);
//        pNew(2) = r*cos(phi);
//        pVec = pNew.normalized();
//        if (std::abs(pVec.dot(normT)) < .999)
//        {
//            // Ensures pnts aren't added very near the wake plane
//            survPnts.row(i) = pNew + trefftzPnt;
//            if (pVec.dot(normal) >= 0)
//            {
//                upperFlag(i) = true;
//            }
//            else
//            {
//                upperFlag(i) = false;
//            }
//            i++;
//        }
    }
}

Eigen::Matrix3d wakeLine::randSurveyPnts(Eigen::MatrixXd &survPnts, Eigen::Matrix<bool,Eigen::Dynamic,1> &upperFlag,Eigen::Vector3d &circDir)
{
    Eigen::Matrix3d POIs;
    
    double xmin,ymin,dx,dy;
    xmin = p1(0);
    ymin = p1(1);
    dx = p2(0)-p1(0);
    dy = p2(1)-p1(1);

    Eigen::Vector3d v = (p2-p1);
    Eigen::Vector3d out = v.cross(normal);
    circDir = v/v.norm();
    
    POIs.row(0) = p1+0.01*out+0.01*normal;
    POIs.row(1) = pMid+0.01*out+0.01*normal;
    POIs.row(2) = p2+0.01*out+0.01*normal;
    
    Eigen::Vector3d newPnt;
    double dz = 0.2;
    
    int i = 0;
    while (i < survPnts.rows())
    {
        newPnt = p1 + v * ((double) rand()/(RAND_MAX)) + out * ((double) rand()/(RAND_MAX));
//        newPnt(2) = p1(2) + (dz * ((double) rand()/(RAND_MAX)) - 0.5*dz);
        newPnt(2) = p1(2) + (dz * ((double) rand()/(RAND_MAX)));

        double cosAngle = (newPnt-p1).dot(normal)/((newPnt-p1).norm()*normal.norm());
//        std::cout << newPnt(0) << "\t" << newPnt(1) << "\t" << newPnt(2) << "\t" << cosAngle << std::endl;
        double eps = 0.00001;
        
        if (std::abs(cosAngle) < eps)
        {
            continue;
        }
        else if (cosAngle > 0)
        {
            upperFlag(i) = true;
            survPnts.row(i) = newPnt;
            i++;
        }
        else
        {
            upperFlag(i) = false;
            survPnts.row(i) = newPnt;
            i++;
        }
    }
    return POIs;
}


void wakeLine::horseshoeW(const Eigen::Vector3d &POI, double &Vy, double &Vz)
{
    double DY1,DY2,DZ1,DZ2,RSQ1,RSQ2;
    double RCORE;
//    if (POI(1) >= p1(1) && POI(1) <= p2(1))
//    {
//        RCORE = 0;
//    }
//    else
//    {
//        RCORE = 0.5*(p2-p1).norm();
//    }
    RCORE = 0.2;
    DY1 = POI(1)-p1(1);
    DY2 = POI(1)-p2(1);
    DZ1 = POI(2)-p1(2);
    DZ2 = POI(2)-p2(2);
//    RSQ1 = DY1*DY1+DZ1*DZ1+RCORE*RCORE;
//    RSQ2 = DY2*DY2+DZ2*DZ2+RCORE*RCORE;
    
//    Vy = M_PI/2*getStrength()*(DZ1/RSQ1-DZ2/RSQ2);
//    Vz = M_PI/2*getStrength()*(-DY1/RSQ1+DY2/RSQ2);
//    Vy = 1/(M_PI*2)*getStrength()*(1/DZ1-1/DZ2);
    Vy = 0;
    Vz = 1/(M_PI*2)*getStrength()*(-1/(DY1-RCORE)+1/(DY2+RCORE));
}


