//
//  surface.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__surface__
#define __CPanel__surface__

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "panel.h"
#include "wakePanel.h"
#include "bodyPanel.h"

class surface
{
    typedef Eigen::Vector3i vertices;
    typedef Eigen::MatrixXd coordinates;
    
    std::vector<panel*> panels;
    short surfID;
    
public:
    surface(const int &surfaceID) : surfID(surfaceID) {}
    
    ~surface()
    {
        for (int i=0; i<panels.size(); i++)
        {
            delete panels[i];
        }
    }
    
    surface(const surface& copy) : surfID(copy.surfID)
    {
        for (int i=0; i<copy.panels.size(); i++)
        {
            panels[i] = new panel(*copy.panels[i]);
        }
    }
    
    void addPanel(const vertices &verts, const coordinates &nodes);

    std::vector<panel*> getPanels() const {return panels;}
    int getID() const {return surfID;}
    
};

#endif /* defined(__CPanel__surface__) */
