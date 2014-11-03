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
    int nNodes = 10;
    std::ofstream fid;
    fid.open("testfile.tri");
    if (fid.is_open())
    {
        // Visual of testfile.tri
        //
        //  6____7,8____9_____10
        //  |\    |\    |\    |
        //  | \ 2 | \ 5 | \ 7 |
        //  |  \  |  \  |  \  |
        //  | 1 \ |3,4\ | 6 \ |
        //  |____\|____\|____\|
        //  1     2     3,4   5
        
        fid << nNodes << "\t" << nTris << "\n";
        
        fid << 0 << "\t" << 0 << "\t" << 0 << "\n"
            << 1 << "\t" << 0 << "\t" << 0 << "\n"
            << 2 << "\t" << 0 << "\t" << 0 << "\n"
            << 2 << "\t" << 0 << "\t" << 0 << "\n"
            << 3 << "\t" << 0 << "\t" << 0 << "\n"
            << 0 << "\t" << 1 << "\t" << 0 << "\n"
            << 1 << "\t" << 1 << "\t" << 0 << "\n"
            << 1 << "\t" << 1 << "\t" << 0 << "\n"
            << 2 << "\t" << 1 << "\t" << 0 << "\n"
            << 3 << "\t" << 1 << "\t" << 0 << "\n";
        
        fid << 1 << "\t" << 2 << "\t" << 6 << "\n"
            << 2 << "\t" << 6 << "\t" << 7 << "\n"
            << 3 << "\t" << 7 << "\t" << 2 << "\n"
            << 3 << "\t" << 7 << "\t" << 2 << "\n"
            << 4 << "\t" << 9 << "\t" << 8 << "\n"
            << 5 << "\t" << 9 << "\t" << 4 << "\n"
            << 5 << "\t" << 10 << "\t" << 9 << "\n";
        
        // Fourth panel is a duplicate of the third panel because wake must stem from two panels on the lifting surface sharing the same edge.
        
        fid << 1 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 10002 << "\n" << 10002 << "\n" << 10002;
        
        fid.close();
        
        geometry testGeom("testfile.tri");
        std::vector<surface*> NLsurfs = testGeom.getNonLiftingSurfs();
        TEST_ASSERT(NLsurfs.size() == 1)
        TEST_ASSERT(NLsurfs[0]->getID() == 1);
        
        TEST_ASSERT(NLsurfs[0]->getPanels().size() == 1);
        
        std::vector<liftingSurf*> Lsurfs = testGeom.getLiftingSurfs();
        
        TEST_ASSERT(Lsurfs[0]->getID() == 2)
        
        TEST_ASSERT(Lsurfs[0]->getPanels().size() == 3);
        
        TEST_ASSERT(Lsurfs[0]->getAllPanels().size() == 6);
        
        wake* w = Lsurfs[0]->getWake();
        
        TEST_ASSERT(w->getPanels()[0]->isTEpanel());
        
        
        

        
    }
    else if (!fid.is_open())
    {
        exit (EXIT_FAILURE);
    }    
}

