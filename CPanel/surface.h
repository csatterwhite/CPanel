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

class geometry;
class bodyStreamline;

class surface
{
protected:
    geometry* geom;
    std::vector<bodyPanel*> panels;
    short surfID;
    bool TEflag; //Surface has sharp trailing edges
    bool LSflag; //Surface is a lifting surface
    std::vector<edge*> trailingEdges;
    Eigen::Vector3d rearStagnationPnt(const Eigen::Vector3d &Vinf, double PG, bodyPanel* &p,double &dist);
    std::vector<edge*> getTrailingEdges();
    
public:
    surface(const int &surfaceID,geometry* geom);
    
//    virtual ~surface();
    
//    surface(const surface& copy) : surfID(copy.surfID)
//    {
//        for (int i=0; i<copy.panels.size(); i++)
//        {
//            panels[i] = new bodyPanel(*copy.panels[i]);
//        }
//    }
    
    virtual void addPanel(bodyPanel* bPan);
    
    void setTEflag() {TEflag = true;}
    void setLSflag() {LSflag = true;}

    std::vector<bodyPanel*> getPanels() const {return panels;}
    int getID() const {return surfID;}
    std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> getStreamlineStartPnts(const Eigen::Vector3d &Vinf,double PG);
    bool sharpTE() {return TEflag;}
    bool isLiftingSurf() {return LSflag;}

    
};

#endif /* defined(__CPanel__surface__) */
