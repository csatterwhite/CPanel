//
//  InfluenceTerms.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/6/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "InfluenceTerms.h"

influenceTerms::influenceTerms(const Eigen::MatrixXd &verts, const Eigen::Vector3d &POI)
{
    // Calculates terms needed to compute source and doublet influence. Ref. Katz and Plotkin
    
    // Verts and POI should be in local reference frame
    
    double rows = verts.rows();
    d.resize(rows);
    m.resize(rows);
    r.resize(rows);
    e.resize(rows);
    h.resize(rows);
    
    for (int i=0; i<verts.rows(); i++)
    {
        Eigen::Vector3d p1;
        Eigen::Vector3d p2;
        if (i!=verts.rows()-1)
        {
            p1 = verts.row(i);
            p2 = verts.row(i+1);
            
        }
        else
        {
            p1 = verts.row(i);
            p2 = verts.row(0);
        }
        d(i) = sqrt(pow(p2(0)-p1(0),2)+pow(p2(1)-p1(1),2));
        m(i) = (p2(1)-p1(1))/(p2(0)-p1(0));
        r(i) = sqrt(pow(POI(0)-p1(0),2)+pow(POI(1)-p1(1),2)+pow(POI(2),2));
        e(i) = pow(POI(0)-p1(0),2)+pow(POI(2),2);
        h(i) = (POI(0)-p1(0))*(POI(1)-p1(1));
    }
}