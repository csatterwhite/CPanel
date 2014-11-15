//
//  main.cpp
//  Octree
//
//  Created by Chris Satterwhite on 4/5/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "octree.h"
#include <vector>

class thing
{
public:
    std::array<double,3> center;
    
    thing()
    {
        center[0] = 0;
        center[1] = 0;
        center[2] = 0;
    }
};

class testOctree : public octree<thing>
{
public:
    std::array<double,3> findRefPoint(const thing* obj)
    {
        return (obj->center);
    }
    
    testOctree() : octree() {}
};

int main()
{
    std::vector<thing*> data;
    thing* obj;
    int nX = 10;
    int nY = 10;
    int nZ = 10;
    for (int i=0; i<nX; i++)
    {
        for (int j=0; j<nY; j++)
        {
            for (int k=0; k<nZ; k++)
            {
                obj = new thing;
                obj->center[0] = i;
                obj->center[1] = j;
                obj->center[2] = k;
                data.push_back(obj);
            }
        }
    }
    testOctree myOctree;
    myOctree.addData(data);
    
}

