//
//  VTUfile.h
//  CPanel
//
//  Created by Chris Satterwhite on 11/25/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__VTUfile__
#define __CPanel__VTUfile__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <typeinfo>


struct cellDataArray
{
    std::string name;
    Eigen::MatrixXd data;
};

struct pntDataArray
{
    std::string name;
    Eigen::MatrixXd data;
};

struct piece
{
    Eigen::MatrixXd pnts;
    Eigen::MatrixXi connectivity;
    std::vector<cellDataArray*> cellData;
    std::vector<pntDataArray*> pntData;
};

class VTUfile
{
    std::vector<piece*> pieces;
    
    void printDoubleArray(std::ofstream &f,std::string name,Eigen::MatrixXd array);
    
    void printIntArray(std::ofstream &f,std::string name,Eigen::MatrixXi array);
    
public:
    VTUfile(std::vector<piece*> pieces) : pieces(pieces) {}
    
    void write(std::string filename);
};

#endif /* defined(__CPanel__VTUfile__) */
