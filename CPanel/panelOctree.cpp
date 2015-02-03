//
//  panelOctree.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "panelOctree.h"

Eigen::Vector3d panelOctree::findRefPoint(const panel &obj)
{
    return obj.getCenter();
}

//panel* panelOctree::getClosestPanel(const Eigen::Vector3d &pnt)
//{
//    // Check panel is within octree
//    if (isInsideOctree(pnt))
//    {
//        std::vector<panel*> pans = findNodeContainingPnt(pnt)->getMembers();
//        
//        panel* closestPan;
//        double dist;
//        double minDist = 1000000;
//        for (int i=0; i<pans.size(); i++)
//        {
//            dist = (pans[i]->getCenter()-pnt).norm();
//            if (dist < minDist)
//            {
//                closestPan = pans[i];
//                minDist = dist;
//            }
//        }
//        if (closestPan->isOnPanel(
//    }
//    
//}