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
    
    cellDataArray(std::string name) : name(name) {}
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
    piece* piece;
    std::string name;
    
    void printDoubleArray(std::ofstream &f,std::string name,Eigen::MatrixXd array);
    
    void printIntArray(std::ofstream &f,std::string name,Eigen::MatrixXi array);
    
    void write();
    
public:
    VTUfile(std::string name, struct piece* piece) : name(name), piece(piece)
    {
        write();
    }
    
    std::string getName() {return name;}
};

#endif /* defined(__CPanel__VTUfile__) */
