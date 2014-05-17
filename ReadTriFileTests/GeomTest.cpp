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
    int nTris = 7;
    int nVerts = 8;
    std::ofstream fid;
    fid.open("testfile.tri");
    if (fid.is_open())
    {
        // Visual of testfile.tri
        //
        //  5_____6_____7_____8
        //  |\    |\    |\    |
        //  | \ 2 | \ 5 | \ 7 |
        //  |  \  |  \  |  \  |
        //  | 1 \ |3,4\ | 6 \ |
        //  |____\|____\|____\|
        //  1     2     3     4
        
        fid << nVerts << "\t" << nTris << "\n";
        
        fid << 0 << "\t" << 0 << "\t" << 0 << "\n"
            << 1 << "\t" << 0 << "\t" << 0 << "\n"
            << 2 << "\t" << 0 << "\t" << 0 << "\n"
            << 3 << "\t" << 0 << "\t" << 0 << "\n"
            << 0 << "\t" << 1 << "\t" << 0 << "\n"
            << 1 << "\t" << 1 << "\t" << 0 << "\n"
            << 2 << "\t" << 1 << "\t" << 0 << "\n"
            << 3 << "\t" << 1 << "\t" << 0 << "\n";
        
        fid << 1 << "\t" << 2 << "\t" << 5 << "\n"
            << 2 << "\t" << 6 << "\t" << 5 << "\n"
            << 3 << "\t" << 6 << "\t" << 2 << "\n"
            << 3 << "\t" << 6 << "\t" << 2 << "\n"
            << 3 << "\t" << 7 << "\t" << 6 << "\n"
            << 4 << "\t" << 7 << "\t" << 3 << "\n"
            << 4 << "\t" << 8 << "\t" << 7 << "\n";
        
        // Fourth panel is a duplicate of the third panel because wake must stem from two panels on the lifting surface sharing the same edge.
        
        fid << 1 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 10002 << "\n" << 10002 << "\n" << 10002;
        
        fid.close();
        
        geometry testGeom;
        testGeom.readGeom("testfile.tri");
        surface* surf = testGeom.getSurfaces()[0];
        TEST_ASSERT_MSG(surf->getID() == 1, "Surface does not have assigned surfID");
        
        TEST_ASSERT_MSG(surf->getPanels().size() == 1, "Incorrect number of panels added to surface");
        
        TEST_ASSERT_MSG(countTEpanels(surf) == 0, "Incorrect TE panel set")
        
        surf = testGeom.getSurfaces()[1];
        
        TEST_ASSERT_MSG(surf->getID() == 2, "Surface does not have assigned surfID")
        
        TEST_ASSERT_MSG(surf->getPanels().size() == 3, "Incorrect number of panels added to surface");
        
        TEST_ASSERT_MSG(countTEpanels(surf) == 2, "Incorrect TE panel set");
        
        surf = testGeom.getSurfaces()[2];
        
        TEST_ASSERT_MSG(surf->getID() == 10002, "Surface does not have assigned surfID")
        
        TEST_ASSERT_MSG(surf->getPanels().size() == 3, "Incorrect number of panels added to surface");
        
        TEST_ASSERT_MSG(countTEpanels(surf) == 0, "Incorrect TE panel set");
        

        
    }
    else if (!fid.is_open())
    {
        exit (EXIT_FAILURE);
    }    
}

int GeomTest::countTEpanels(surface *surf)
{
    int count = 0;
    std::vector<panel*> panels = surf->getPanels();
    for (int i=0; i<panels.size(); i++)
    {
        if (panels[i]->isTEpanel())
        {
            count++;
        }
    }
    return count;
}