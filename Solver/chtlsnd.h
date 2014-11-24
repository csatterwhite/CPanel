//
//  chtlsnd.h
//  CPanel
//
//  Created by Chris Satterwhite on 11/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__chtlsnd__
#define __CPanel__chtlsnd__

#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/QR>

class chtlsnd
{
    Eigen::MatrixXd F;
    
    
public:
    chtlsnd(const Eigen::Vector3d &X0, const Eigen::MatrixXd &Xf, int order, const Eigen::MatrixXd &Xb, const Eigen::MatrixXd &Vb, Eigen::Vector3d V0);
    
};

#endif /* defined(__CPanel__chtlsnd__) */
