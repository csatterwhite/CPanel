//
//  geometryTests.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/22/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__geometryTests__
#define __CPanel__geometryTests__

#include <stdio.h>
#include <iostream>
#include "cpptest/cpptest.h"
#include "Eigen/Dense"
#include <fstream>
#include "geometry.h"
#include "surface.h"
#include "liftingSurf.h"
#include "wake.h"
#include "bodyPanel.h"
#include "wakePanel.h"

class GeomTests : public Test::Suite
{
    std::string testTriFile;
    geometry* testGeom;
    void test_readTri();
    void test_neighborSearch();
    void test_createOctree();
    void makeTestTriFile();
    int nTris,nNodes;
    
    // Visual of testfile.tri
    // Representative of wing body junction with a shed wake.  Addresses all special cases in neighbor search.  Panel numbers prefixed with "B" are body panels and "W" are wake panels
    //
    // B2, B5, and B10 are lower panels on opposite side of wake
    //
    // 0   1,2____3,4_____5______6
    // .    |\     |\     |     /|
    // .    | \ B3 | \ B6 | B7 / |
    // .    |  \   |  \   |   /  |
    // .    | B1\  | B4\  |  /   |
    // .    | B2 \ | B5 \ | / B8 |
    // 1  7,14____8,15___9,16____10,11
    // .    |\     |\     |     /|
    // .    | \ W13| \W15 | B9 / |
    // .    |  \   |  \   |B10/  |
    // .    |   \  |   \  |  /   |
    // .    |W12 \ |W14 \ | / B11|
    // 2    |_____\|_____\|/_____|
    // |    17    18     12,19   13
    // x
    //   y- -3....-2.....-1......0
    
public:
    GeomTests();
    
    ~GeomTests()
    {
        delete testGeom;
    }
};

#endif /* defined(__CPanel__geometryTests__) */
