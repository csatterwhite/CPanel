//
//  main.cpp
//  Solver
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "geometry.h"
#include "runCase.h"

int main(int argc, const char * argv[])
{
    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
    std::string inFile = "genericAC.tri";
    std::string outFile = "genericAC_neighborCheck_e3.vtu";

    time_t ts;
    time(&ts);
    time_t tf;
    geometry temp(path+inFile);
    
    std::cout <<  temp.getOctree()->getMembers().size() << std::endl;
    
//    std::string neighborFile = "NeighborCheck_genericAC.txt";
//    std::ofstream fid;
//    fid.open(path+neighborFile);
//    std::vector<panel*> pans = temp.getPanels();
//    Eigen::MatrixXd nodes = temp.getNodes();
//    Eigen::MatrixXi conn;
//    Eigen::VectorXi nNeighbs;
//    for (int i=0; i<pans.size(); i++)
//    {
//        conn.row(i) = pans[i]->getVerts();
//        nNeighbs(i) = pans[i]->getNeighbors().size();
//    }
//    fid << nodes.rows() << "\t" << conn.rows() << "\n";
//    for (int i=0; i<nodes.rows(); i++)
//    {
//        fid << nodes(i,0) << "\t" << nodes(i,1) << "\t" << nodes(i,2) << "\n";
//    }
//    for (int i=0; i<conn.rows(); i++)
//    {
//        fid << conn(i,0) << "\t" << conn(i,1) << "\t" << conn(i,2) << "\n";
//    }
//    for (int i=0; i<conn.rows(); i++)
//    {
//        fid << nNeighbs(i,0) << "\n";
//    }
//    fid.close();
    
    runCase case1(&temp,1,0,0,path+outFile);
    
    time(&tf);
    std::cout << "Elapsed time for program execution : " << difftime(tf,ts) << " seconds" << std::endl;
}
