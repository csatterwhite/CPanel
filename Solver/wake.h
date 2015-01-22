//
//  wake.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wake__
#define __CPanel__wake__

#include <stdio.h>
#include <vector>
#include <cmath>
#include "wakePanel.h"
#include "wakeLine.h"

class wakePanel;

class wake
{
    std::vector<wakePanel*> wpanels;
    std::vector<wakeLine*> wakeLines;
    std::vector<wakePanel*> vortexSheets;
    Eigen::MatrixXd* nodes;
    double yMin;
    double yMax;
    double x0,xf,z0,zf;
    
    
    void setWakeDimensions();
    short edgeVerts(wakePanel* p);
    wakeLine* findWakeLine(double y);
    double Vradial(Eigen::Vector3d pWake);
    Eigen::Vector3d pntInWake(double x, double y);
    Eigen::Vector3d pntVel(Eigen::Matrix<double,1,3> POI, Eigen::MatrixXd pntCloud,Eigen::Matrix<bool,Eigen::Dynamic,1> upperFlag, Eigen::Vector3d Vinf);
    
public:
    wake(Eigen::MatrixXd* nodes) : nodes(nodes), yMin(0) {}
    
    ~wake();
    
    wake(const wake& copy);
    
    void addPanel(wakePanel* wPan);
    
    std::vector<wakePanel*> getPanels() const {return wpanels;}
    
    void setNeighbors(panelOctree* oct);
    
    void trefftzPlane(double Vinf,double Sref, double &CL, double &CD, Eigen::VectorXd &yLoc, Eigen::VectorXd &Cl, Eigen::VectorXd &Cd);
    
    std::vector<wakeLine*> getWakeLines() {return wakeLines;}
    
    double wakeStrength(double y);
    
    std::vector<wakePanel*> getVortexSheets() {return vortexSheets;}
};

#endif /* defined(__CPanel__wake__) */
