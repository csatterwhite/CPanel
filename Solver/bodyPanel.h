//
//  bodyPanel.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__bodyPanel__
#define __CPanel__bodyPanel__

#include <iostream>
#include "panel.h"

class bodyPanel : public panel
{
    double sigma;
    double mu;
    bool TEpanel;
    
public:
    bodyPanel(const Eigen::VectorXi &panelVertices,Eigen::MatrixXd* nodes,int surfID) : panel(panelVertices,nodes,surfID), TEpanel(false)
    {
        sigma = 1;
        mu = 1;
    }
    
    
    bodyPanel(const bodyPanel &copy) : panel(copy.verts,copy.nodes,copy.ID), sigma(copy.sigma), mu(copy.mu), TEpanel(copy.TEpanel)
    {
    }
    
    void setSigma(Eigen::Vector3d Vinf, double Vnorm)
    {
        sigma = -Vinf.dot(normal)+Vnorm;
    }
    
    void setMu(double doubletStrength)
    {
        mu = doubletStrength;
    }
    
    void setPotential(Eigen::Vector3d Vinf)
    {
        potential = mu+Vinf.dot(center);
    }
    
    void setTEpanel() {TEpanel = true;}
    
    double getSigma() {return sigma;}
    double getMu() {return mu;}
    bool isTEpanel() {return TEpanel;}
};

#endif /* defined(__CPanel__bodyPanel__) */
