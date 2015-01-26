//
//  geometryTests.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/22/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "geometryTests.h"

GeomTests::GeomTests()
{
    makeTestTriFile();
    TEST_ADD(GeomTests::test_readTri);
    TEST_ADD(GeomTests::test_neighborSearch);
    TEST_ADD(GeomTests::test_createOctree);
}

void GeomTests::test_readTri()
{
    // Tests createSurvaces method in geometry class, as well as addPanel method in the surface class.
    // Create Geometry
    Eigen::Vector3d cg = Eigen::Vector3d::Zero();
    testGeom = new geometry(testTriFile,false,1,1,1,cg);
    
    // 5 body panels should be created in a nonlifting surface with surfID of 2;
    std::vector<surface*> NLsurfs = testGeom->getNonLiftingSurfs();
    TEST_ASSERT(NLsurfs.size() == 1);
    TEST_ASSERT(NLsurfs[0]->getID() == 2);
    TEST_ASSERT(NLsurfs[0]->getPanels().size() == 5);
    
    // 6 body panels should be created in a lifting surface with surfID of 1.  4 corresponding wake panels should have been added to its wake
    std::vector<liftingSurf*> Lsurfs = testGeom->getLiftingSurfs();
    TEST_ASSERT(Lsurfs.size() == 1);
    TEST_ASSERT(Lsurfs[0]->getID() == 1);
    TEST_ASSERT(Lsurfs[0]->getPanels().size() == 6); // Just Wing Panels
    TEST_ASSERT(Lsurfs[0]->getWakePanels().size() == 4);
    TEST_ASSERT(Lsurfs[0]->getAllPanels().size() == 10); // Wing+Wake Panels
}

void GeomTests::test_neighborSearch()
{
    std::vector<bodyPanel*> *bPanels = testGeom->getBodyPanels();
    Eigen::VectorXi nNeighbs(11);
    nNeighbs << 2,1,2,3,1,2,2,2,3,1,1;
    for (int i=0; i<(*bPanels).size(); i++)
    {
        TEST_ASSERT_MSG((*bPanels)[i]->getNeighbors().size() == nNeighbs(i), "Incorrect Number of Neighboring Panels Set" );
    }
    
    std::vector<wakePanel*> *wPanels = testGeom->getWakePanels();
    Eigen::VectorXi TEpans(4);
    TEpans << 0,1,0,1;
    for (int i=0; i<wPanels->size(); i++)
    {
        TEST_ASSERT_MSG((*wPanels)[i]->isTEpanel() == TEpans(i), "Trailing Edge Panels Not Set Properly");
    }
    
    TEST_ASSERT_MSG((*wPanels)[1]->getUpper() == (*bPanels)[0], "Upper Parent Panel not set properly")
    
    TEST_ASSERT_MSG((*wPanels)[1]->getLower() == (*bPanels)[1], "Lower Parent Panel not set properly")
    
    TEST_ASSERT_MSG((*wPanels)[3]->getUpper() == (*bPanels)[3], "Upper Parent Panel not set properly at wing-body joint")
    
    TEST_ASSERT_MSG((*wPanels)[3]->getLower() == (*bPanels)[4], "Lower Parent Panel not set properly at wing-body joint")
}

void GeomTests::test_createOctree()
{
    TEST_ASSERT(testGeom->getOctree()->getMembers().size() == nTris);
}

void GeomTests::makeTestTriFile()
{
    // Construct tri file to test
    std::ofstream fid;
    testTriFile = "testfile.tri";
    fid.open(testTriFile);
    if (fid)
    {
        // See Visual of testfile.tri in geometryTests.h
        // Nodes with z=-1 (2,4,11) are for lower panel tris.
        
        nTris = 15;
        nNodes = 19;
        double eps = pow(10,-7); // Addresses floating point error in wake nodes intended to be coincident with body nodes seen in .tri files output by OpenVSP.
        
        fid << nNodes << "\t" << nTris << "\n";
        
        fid << 0 << "\t" << -3 << "\t" << 0 << "\n"
        << 0 << "\t" << -3 << "\t" << -1 << "\n"
        << 0 << "\t" << -2 << "\t" << 0 << "\n"
        << 0 << "\t" << -2 << "\t" << -1 << "\n"
        << 0 << "\t" << -1 << "\t" << 0 << "\n"
        << 0 << "\t" << 0 << "\t" << 0 << "\n"
        << 1 << "\t" << -3 << "\t" << 0 << "\n"
        << 1 << "\t" << -2 << "\t" << 0 << "\n"
        << 1 << "\t" << -1 << "\t" << 0 << "\n"
        << 1 << "\t" << 0 << "\t" << 0 << "\n"
        << 1 << "\t" << 0 << "\t" << -1 << "\n"
        << 2 << "\t" << -1 << "\t" << 0 << "\n"
        << 2 << "\t" << 0 << "\t" << 0 << "\n"
        << 1 << "\t" << -3 << "\t" << eps << "\n"
        << 1 << "\t" << -2 << "\t" << eps << "\n"
        << 1 << "\t" << -1 << "\t" << eps << "\n"
        << 2 << "\t" << -3 << "\t" << 0 << "\n"
        << 2 << "\t" << -2 << "\t" << 0 << "\n"
        << 2 << "\t" << -1 << "\t" << 0 << "\n";
        
        fid << 1 << "\t" << 7 << "\t" << 8 << "\n"
        << 2 << "\t" << 8 << "\t" << 7 << "\n"
        << 1 << "\t" << 8 << "\t" << 3 << "\n"
        << 3 << "\t" << 8 << "\t" << 9 << "\n"
        << 4 << "\t" << 9 << "\t" << 8 << "\n"
        << 3 << "\t" << 9 << "\t" << 5 << "\n"
        << 5 << "\t" << 9 << "\t" << 6 << "\n"
        << 6 << "\t" << 9 << "\t" << 10 << "\n"
        << 10 << "\t" << 9 << "\t" << 12 << "\n"
        << 11 << "\t" << 12 << "\t" << 9 << "\n"
        << 10 << "\t" << 12 << "\t" << 13 << "\n"
        << 14 << "\t" << 17 << "\t" << 18 << "\n"
        << 14 << "\t" << 18 << "\t" << 15 << "\n"
        << 15 << "\t" << 18 << "\t" << 19 << "\n"
        << 15 << "\t" << 19 << "\t" << 16 << "\n";
        
        fid << 1 << "\n" << 1 << "\n" << 1 << "\n" << 1 << "\n" << 1 << "\n" << 1 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 10001 << "\n" << 10001 << "\n" << 10001 << "\n" << 10001;
        
        fid.close();
    }
}