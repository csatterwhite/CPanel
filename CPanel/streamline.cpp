//
//  streamline.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/26/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "streamline.h"

streamline::streamline(const Eigen::Vector3d &startPnt, double xMax, double tol, const Eigen::Vector3d &Vinf, double PG, geometry* geom) : Vinf(Vinf), PG(PG), geom(geom)
{
    coeff5 << 16.0/135 , 0 , 6656.0/12825 , 28561.0/56430 , -9.0/50 , 2.0/55;
    coeff4 << 25.0/216 , 0 , 1408.0/2565 , 2197.0/4104 , -0.2 , 0;
    
    double error;
    double dt = 0.05/Vinf.norm();
    Eigen::Vector3d nextPnt;
    nextPnt = rkf(startPnt, dt, error);
    
    while (nextPnt(0) < xMax && pnts.size() < 600)
    {
        if (error == 0)
        {
            dt = dt;
        }
        else
        {
            dt = dt*pow((tol*dt/(2*error)),0.25);
        }
        nextPnt = rkf(nextPnt, dt, error);
        if (nextPnt(0) > xMax)
        {
            int dummy = 1;
        }
    }
}

Eigen::Vector3d streamline::rkf(const Eigen::Vector3d &x0,double dt,double &error)
{
    Eigen::Vector3d step4,step5;
    Eigen::Matrix<double,6,3> ks;
    Eigen::Vector3d x2,x3,x4,x5,x6;
    Eigen::Matrix<double,1,3> x1;
    x1(0) = x0(0);
    x1(1) = x0(1);
    x1(2) = x0(2);
    
    Eigen::Vector3d vPnt = geom->pntVelocity(x0,Vinf,PG);
    velocities.push_back(vPnt);
    pnts.push_back(x0);
    
    ks.row(0) = dt*vPnt;
    
    x2 = x1 + 1.0/4*ks.row(0);
    ks.row(1) = dt*geom->pntVelocity(x2,Vinf,PG);
    
    x3 = x1 + 1.0/32 * (3*ks.row(0) + 9.0*ks.row(1));
    ks.row(2) = dt*geom->pntVelocity(x3,Vinf,PG);
    
    x4 = x1 + 1.0/2197 * (1932*ks.row(0) - 7200*ks.row(1) + 7296*ks.row(2));
    ks.row(3) = dt*geom->pntVelocity(x4,Vinf,PG);
    
    x5 = x1 + 439.0/216*ks.row(0) - 8*ks.row(1) + 3680.0/513*ks.row(2) - 845.0/4104*ks.row(3);
    ks.row(4) = dt*geom->pntVelocity(x5,Vinf,PG);
    
    x6 = x1 - 8.0/27*ks.row(0) + 2*ks.row(1) - 3544.0/2565*ks.row(2) +1859.0/4104*ks.row(3) - 11.0/40*ks.row(4);
    ks.row(5) = dt*geom->pntVelocity(x6,Vinf,PG);
    
    step5 = coeff5*ks;
    step4 = coeff4*ks;
    error = (step5-step4).norm();
    
    return x0+step5;
}