//
//  InfluenceTerms.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/6/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__InfluenceTerms__
#define __CPanel__InfluenceTerms__

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>

struct influenceTerms
{
    Eigen::VectorXd d,m,r,e,h;
    influenceTerms(const Eigen::MatrixXd &verts, const Eigen::Vector3d &POI);

};

#endif /* defined(__CPanel__InfluenceTerms__) */
