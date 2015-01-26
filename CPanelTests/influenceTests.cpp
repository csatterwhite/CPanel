//
//  influenceTests.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/24/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "influenceTests.h"

influenceTests::influenceTests()
{
    createTestPanel();
        TEST_ADD(influenceTests::testPntClose)
        TEST_ADD(influenceTests::testPntFar)
//        TEST_ADD(influenceTests::testCollocationPnt)
    //    TEST_ADD(influenceTests::testDoublet)
    //    TEST_ADD(influenceTests::testSource)
}

void influenceTests::createTestPanel()
{
    testNodes.resize(3,3);
    testNodes(0,0) = 1;
    testNodes(0,1) = 0;
    testNodes(0,2) = 0;
    testNodes(1,0) = -1;
    testNodes(1,1) = 1;
    testNodes(1,2) = 0;
    testNodes(2,0) = -1;
    testNodes(2,1) = -1;
    testNodes(2,2) = 0;
    Eigen::Vector3i verts;
    verts << 0,1,2;
    Eigen::Vector3d bezNorm;
    bezNorm << 0,0,1;
    std::vector<edge*> edges;
    testPan = new bodyPanel(verts,&testNodes,edges,bezNorm,1,false);
}

void influenceTests::testPntClose()
{
    Eigen::Vector3d POI;
    POI << 2,2,2; // |POI|/(side length) ~ 1.2 so panel formulation will be used;

    double phiSrc,phiDub;
    Eigen::Vector3d vSrc,vDub;
    phiSrc = 0;
    phiDub = 0;
    vSrc = Eigen::Vector3d::Zero();
    vDub = Eigen::Vector3d::Zero();
    testPan->panelPhiInf(POI, phiSrc, phiDub);
    testPan->panelVInf(POI, vSrc, vDub);
    TEST_ASSERT(isEqual(phiSrc,0.043411,5))
    TEST_ASSERT(isEqual(vSrc(0),0.007363,5))
    TEST_ASSERT(isEqual(vSrc(1),0.006451,5))
    TEST_ASSERT(isEqual(vSrc(2),0.006645,5))

    TEST_ASSERT(isEqual(phiDub,0.006645,5))
    TEST_ASSERT(isEqual(vDub(0),-0.003381,5))
    TEST_ASSERT(isEqual(vDub(1),-0.003026,5))
    TEST_ASSERT(isEqual(vDub(2),0.000158,5))
}

void influenceTests::testPntFar()
{
    Eigen::Vector3d POI;
    POI << 7,7,7; // |POI|/(side length) ~ 5.5 so far field approximation will be used;

    double phiSrc,phiDub;
    Eigen::Vector3d vSrc,vDub;
    phiSrc = 0;
    phiDub = 0;
    vSrc = Eigen::Vector3d::Zero();
    vDub = Eigen::Vector3d::Zero();
    testPan->panelPhiInf(POI, phiSrc, phiDub);
    testPan->panelVInf(POI, vSrc, vDub);
    TEST_ASSERT(isEqual(phiSrc,0.012919,5))
    TEST_ASSERT(isEqual(vSrc(0),0.000625,5))
    TEST_ASSERT(isEqual(vSrc(1),0.000596,5))
    TEST_ASSERT(isEqual(vSrc(2),0.000596,5))
    
    TEST_ASSERT(isEqual(phiDub,0.0005958,5))
    TEST_ASSERT(isEqual(vDub(0),-0.000086,6))
    TEST_ASSERT(isEqual(vDub(1),-0.000082,6))
    TEST_ASSERT(isEqual(vDub(2),0.0000027,6))
}

//void influenceTests::testCollocationPnt()
//{
//    Eigen::Vector3d POI;
//    POI << -1.0/3,0,0;
//
//    double phiSource;
//    Eigen::Vector3d vSource;
//    phiSource = p.sourcePhi(1,POI,vertsLocal);
//    vSource = p.sourceV(1,POI,vertsLocal);
//    TEST_ASSERT(phiSource == 0)
//    TEST_ASSERT(vSource(0) == 0)
//    TEST_ASSERT(vSource(1) == 0)
//    TEST_ASSERT(vSource(2) == 0.5)
//
//    double phiDoublet;
//    Eigen::Vector3d vDoublet;
//    phiDoublet = p.doubletPhi(1,POI,vertsLocal);
//    vDoublet = p.doubletV(1,POI,vertsLocal);
//    TEST_ASSERT(phiDoublet == -0.5)
//    TEST_ASSERT(vDoublet(0) == 0)
//    TEST_ASSERT(vDoublet(1) == 0)
//    TEST_ASSERT(isEqual(vDoublet(2),1.432394487,9))
//}

bool influenceTests::isEqual(double var1, double var2, int decimalPrecision)
{
    double eps = pow(10,-decimalPrecision);
    return (std::abs(var1-var2)<eps);
    
}

//void influenceTests::testDoublet()
//{
//    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
//    std::string filename = "doubletCheck.txt";
//    std::ifstream fid;
//    fid.open(path+filename);
//    if (fid.is_open())
//    {
//        Eigen::MatrixXd nodes;
//        nodes.resize(3,3);
//        Eigen::Vector3i verts;
//        verts << 0,1,2;
//        Eigen::Vector3d col;
//        Eigen::MatrixXd pnts(9,3);
//        Eigen::VectorXd pot(9);
//        Eigen::MatrixXd vel(9,3);
//
//        for (int i=0; i<3; i++)
//        {
//            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
//        }
//        fid >> col(0) >> col(1) >> col(2);
//        for (int i=0; i<9; i++)
//        {
//            fid >> pnts(i,0) >> pnts(i,1) >> pnts(i,2);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> pot(i);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> vel(i,0) >> vel(i,1) >> vel(i,2);
//        }
//        vel *= -1; //Different sign convention
//        Eigen::VectorXd phiCalc(9);
//        Eigen::MatrixXd Vcalc(9,3);
//        bodyPanel pan(verts,&nodes,1);
//        Eigen::Matrix3d localVerts = pan.getLocalVerts();
//
//        for (int i=0; i<9; i++)
//        {
//            phiCalc(i) = pan.doubletPhi(1,pnts.row(i),localVerts);
//            TEST_ASSERT(isEqual(phiCalc(i),pot(i),5))
//            Vcalc.row(i) = pan.doubletV(true, pnts.row(i));
//            for (int j=0; j<3; j++)
//            {
//                std::cout << Vcalc(i,j) << "," << vel(i,j) << std::endl;
//                TEST_ASSERT(isEqual(Vcalc(i,j),vel(i,j),5))
//            }
//        }
//    }
//    else
//    {
//        std::cout << "ERROR : File not found" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//}
//
//void influenceTests::testSource()
//{
//    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
//    std::string filename = "sourceCheck.txt";
//    std::ifstream fid;
//    fid.open(path+filename);
//    if (fid.is_open())
//    {
//        Eigen::MatrixXd nodes;
//        nodes.resize(3,3);
//        Eigen::Vector3i verts;
//        verts << 0,1,2;
//        Eigen::Vector3d col;
//        Eigen::MatrixXd pnts(9,3);
//        Eigen::VectorXd pot(9);
//        Eigen::MatrixXd vel(9,3);
//
//        for (int i=0; i<3; i++)
//        {
//            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
//        }
//        fid >> col(0) >> col(1) >> col(2);
//        for (int i=0; i<9; i++)
//        {
//            fid >> pnts(i,0) >> pnts(i,1) >> pnts(i,2);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> pot(i);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> vel(i,0) >> vel(i,1) >> vel(i,2);
//        }
//        pot *= -1; //Different sign convention;
//        Eigen::VectorXd phiCalc(9);
//        Eigen::MatrixXd Vcalc(9,3);
//        bodyPanel pan(verts,&nodes,1);
//        Eigen::Matrix3d localVerts = pan.getLocalVerts();
//
//        for (int i=0; i<9; i++)
//        {
//            phiCalc(i) = pan.sourcePhi(1,pnts.row(i),localVerts);
//            TEST_ASSERT(isEqual(phiCalc(i),pot(i),5))
//            Vcalc.row(i) = pan.sourceV(true, pnts.row(i));
//            for (int j=0; j<3; j++)
//            {
//                TEST_ASSERT(isEqual(-Vcalc(i,j),vel(i,j),5))
//            }
//        }
//    }
//    else
//    {
//        std::cout << "ERROR : File not found" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//}
