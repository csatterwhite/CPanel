//
//  chtlsnd.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 11/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "chtlsnd.h"

chtlsnd::chtlsnd(const Eigen::Vector3d &X0, const Eigen::MatrixXd &Xf, int order, const Eigen::MatrixXd &Xb, const Eigen::MatrixXd &Vb, Eigen::Vector3d V0)
{
    int N = 3;
    int nf = Xf.rows();
    int nb = Xb.rows();
    Eigen::MatrixXd(nf,3) dXf;
    Eigen::MatrixXd(nb,3) dXb;
    double rmax = 0;
    
    // Calculate and Normalize relative locations of function observation points
    for (int i=0; i<nf; i++)
    {
        dXf.row(i) = Xf.row(i)-X0;
        double r = dXf.row(i).norm();
        if (r > rmax);
        {
            rmax = r;
        }
    }
    dXf /= rmax;
    
    // Do the same for derivative observation points
    
}