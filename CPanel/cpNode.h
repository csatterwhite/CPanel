//
//  cpNode.h
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 2/2/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__cpNode__
#define __CPanel___Unstructured_Panel_Code__cpNode__

#include <stdio.h>
#include <Eigen/Dense>
#include <vector>

class edge;
class bodyPanel;

class cpNode
{
    Eigen::Vector3d pnt;
    int index;
    std::vector<edge*> edges;
    std::vector<bodyPanel*> bodyPans;
    
    bool TEnode;
    
public:
    cpNode(Eigen::Vector3d pnt,int index);
    
//    cpNode(const cpNode& copy);
    
    Eigen::Vector3d operator-=(const cpNode &rhs);
    
    Eigen::Vector3d operator+=(const cpNode &rhs);
    
    void addEdge(edge* e);
    void addBodyPanel(bodyPanel* p);
    
    edge* getTE(edge* exception);
    
    void setTE();
    void setIndex(int i);
    
    Eigen::Vector3d getPnt() const {return pnt;}
    
    int getIndex() const {return index;}
    
    std::vector<edge*> getEdges() {return edges;}
    
    std::vector<bodyPanel*> getBodyPans() {return bodyPans;}
    
    edge* getOtherTrailEdge(edge* current);
    
    bool isTE() {return TEnode;}
    
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__cpNode__) */
