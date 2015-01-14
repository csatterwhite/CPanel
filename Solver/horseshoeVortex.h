//
//  horseshoeVortex.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/10/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__horseshoeVortex__
#define __CPanel__horseshoeVortex__

#include <stdio.h>
#include <Eigen/Dense>
#include "bodyPanel.h"



class horseshoeVortex
{
    double gamma;
    Eigen::Vector3d p1,p2;
    bodyPanel *upper,*lower;
    
public:
    horseshoeVortex(Eigen::Vector3d p1, Eigen::Vector3d p2, bodyPanel* upper, bodyPanel* lower) : p1(p1), p2(p2), upper(upper), lower(lower) {}
    
    double phiInfluence(Eigen::Vector3d POI);
    double horseshoePhi(Eigen::Vector3d POI);
    void setStrength() {gamma = upper->getMu()-lower->getMu();}
    double getStrength() {return gamma;}
    bodyPanel* getUpper() {return upper;}
    bodyPanel* getLower() {return lower;}
    
    
};
#endif /* defined(__CPanel__horseshoeVortex__) */
