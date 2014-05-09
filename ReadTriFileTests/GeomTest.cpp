//
//  GeomTest.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 3/27/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "GeomTest.h"

GeomTest::GeomTest()
{
    TEST_ADD(GeomTest::test_readGeom);
}

void GeomTest::test_readGeom()
{
   // Construct tri file to test
    int nTris = 6;
    int nVerts = 8;
    std::ofstream fid;
    fid.open("testfile.tri");
    if (fid.is_open())
    {
        fid << nVerts << "\t" << nTris << "\n";
        
        fid << 0 << "\t" << 0 << "\t" << 0 << "\n"
            << 1 << "\t" << 0 << "\t" << 0 << "\n"
            << 2 << "\t" << 0 << "\t" << 0 << "\n"
            << 3 << "\t" << 0 << "\t" << 0 << "\n"
            << 0 << "\t" << 1 << "\t" << 0 << "\n"
            << 1 << "\t" << 1 << "\t" << 0 << "\n"
            << 2 << "\t" << 1 << "\t" << 0 << "\n"
            << 3 << "\t" << 1 << "\t" << 0 << "\n";
        
        fid << 1 << "\t" << 5 << "\t" << 2 << "\n"
            << 2 << "\t" << 5 << "\t" << 6 << "\n"
            << 3 << "\t" << 2 << "\t" << 6 << "\n"
            << 3 << "\t" << 6 << "\t" << 7 << "\n"
            << 4 << "\t" << 3 << "\t" << 7 << "\n"
            << 4 << "\t" << 7 << "\t" << 8 << "\n";
        
        fid << 1 << "\n" << 2 << "\n" << 2 << "\n" << 10001 << "\n" << 10001 << "\n" << 10001;
        
        fid.close();
        
        geometry testGeom;
        testGeom.readGeom("testfile.tri");
        TEST_ASSERT_MSG(testGeom.getSurfaces()[0]->getID() == 1, "Surface does not have assigned surfID");
        
        TEST_ASSERT_MSG(testGeom.getSurfaces()[1]->getID() == 2, "Surface does not have assigned surfID");
        
        TEST_ASSERT_MSG(testGeom.getSurfaces()[1]->getID() == 2, "Surface does not have assigned surfID");
        
        TEST_ASSERT_MSG(testGeom.getSurfaces()[2]->getID() == 10001, "Surface does not have assigned surfID");
        
        TEST_ASSERT_MSG(testGeom.getSurfaces()[2]->getID() == 10001, "Surface does not have assigned surfID");
        
        TEST_ASSERT_MSG(testGeom.getSurfaces()[2]->getID() == 10001, "Surface does not have assigned surfID");
        
    }
    else if (!fid.is_open())
    {
        exit (EXIT_FAILURE);
    }    
}