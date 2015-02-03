//
//  wakeLine.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wakeLine__
#define __CPanel__wakeLine__

#include <stdio.h>
#include <cmath>
#include "bodyPanel.h"
#include "edge.h"
#include "cpNode.h"

class bodyPanel;

class wakeLine
{
    bodyPanel* upper;
    bodyPanel* lower;
    Eigen::Vector3d p1;
    Eigen::Vector3d pMid;
    Eigen::Vector3d p2;
    Eigen::Vector3d normal;
    
    void setDimensions();
    
public:
    wakeLine(bodyPanel* upper, bodyPanel* lower,Eigen::Vector3d normal);
    
    void surveyPnts(double y,double distDownStream, double x0Wake, double z0Wake, Eigen::Vector3d &trefftzPnt, Eigen::MatrixXd &survPnts, Eigen::Matrix<bool,Eigen::Dynamic,1> &upperFlag);
    Eigen::Matrix3d randSurveyPnts(Eigen::MatrixXd &survPnts, Eigen::Matrix<bool,Eigen::Dynamic,1> &upperFlag, Eigen::Vector3d &circDir);
    
    double getY() {return pMid(1);}
    bodyPanel* getUpper() {return upper;}
    bodyPanel* getLower() {return lower;}
    Eigen::Vector3d getP1() {return p1;}
    Eigen::Vector3d getPMid() {return pMid;}
    Eigen::Vector3d getP2() {return p2;}
    Eigen::Vector3d getNormal() {return normal;}
    double getStrength();
    void horseshoeW(const Eigen::Vector3d &POI, double &Vy, double &Vz);
};

#endif /* defined(__CPanel__wakeLine__) */
