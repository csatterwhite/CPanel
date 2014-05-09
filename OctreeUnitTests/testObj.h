//
//  testObj.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/8/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__testObj__
#define __CPanel__testObj__

#include <iostream>
#include <array>

class testObj
{
    std::array<double,3> center;
public:
    testObj(std::array<double,3> point) : center(point) {}
    
    std::array<double,3> getCenter() const {return center;}
    
};

#endif /* defined(__CPanel__testObj__) */
