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
    Eigen::MatrixXd G;
    Eigen::MatrixXd H;
    
    Eigen::MatrixXi derivSequence(int q, int N);
    Eigen::MatrixXi sortBySum(Eigen::MatrixXi m);
    Eigen::MatrixXi insertRow(const Eigen::Matrix<int,Eigen::Dynamic,3> &m, const Eigen::Matrix<int,1,3> &insert, int row);
    int factorial(int i);
    
public:
    chtlsnd(const Eigen::Matrix<double,1,3> &X0, const Eigen::Matrix<double,Eigen::Dynamic,3> &Xf, int order, const Eigen::Matrix<double,Eigen::Dynamic,3> &Xb, const Eigen::Matrix<double,Eigen::Dynamic,3> &Vb, Eigen::Vector3d V0);
    
};

#endif /* defined(__CPanel__chtlsnd__) */
