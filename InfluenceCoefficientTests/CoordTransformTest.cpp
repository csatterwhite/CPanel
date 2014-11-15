//
//  CoordTransformTest.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 5/19/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "CoordTransformTest.h"

//CoordTransformTest::CoordTransformTest()
//{
//    TEST_ADD(CoordTransformTest::testTransform);
//}
//
//void CoordTransformTest::testTransform()
//{
//    Eigen::Vector3d vec = {1,1,1};
//    Eigen::Matrix3d mat1;
//    Eigen::Matrix3d mat2;
//    Eigen::Matrix3d mat3;
//    Eigen::Matrix3d mat4;
//    
//    mat1 << 1,0,0, 0,1,0, 0,0,1;
//    mat2 << 0,1,0, -1,0,0, 0,0,1; //90 degree rotation about Z
//    mat3 << 0,0,1, 0,1,0, -1,0,0; //90 degree rotation about Y
//    mat4 << 1,0,0, 0,0,1, 0,-1,0; //90 degree rotation about X
//    
//    panel p(1);
//    Eigen::Vector3d newVec = p.transformCoordinates(vec, mat1, mat2);
//    
//    TEST_ASSERT(newVec(0) == 1 && newVec(1) == -1 && newVec(2) == 1);
//    
//    newVec = p.transformCoordinates(vec, mat1, mat3);
//    
//    TEST_ASSERT(newVec(0) == 1 && newVec(1) == 1 && newVec(2) == -1);
//    
//    newVec = p.transformCoordinates(vec, mat1, mat4);
//    
//    TEST_ASSERT(newVec(0) = 1 && newVec(1) == 1 && newVec(2) == -1);
//    
//}