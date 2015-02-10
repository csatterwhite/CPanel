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
    
    pntDataArray(std::string name) : name(name) {}
};

struct piece
{
    Eigen::MatrixXd pnts;
    Eigen::MatrixXi connectivity;
    std::vector<cellDataArray> cellData;
    std::vector<pntDataArray> pntData;
};

class VTUfile
{
    std::vector<piece> pieces;
    std::string name;
    
    void printDoubleArray(std::ofstream &f,std::string name,Eigen::MatrixXd array);
    
    void printIntArray(std::ofstream &f,std::string name,Eigen::MatrixXi array);
    
    void write();
    
public:
    VTUfile(std::string name, piece p) : name(name)
    {
        pieces.push_back(p);
        write();
    }
    
    VTUfile(std::string name, std::vector<piece> pieces) : name(name),pieces(pieces)
    {
        write();
    }
    
    std::string getName() {return name;}
};

#endif /* defined(__CPanel__VTUfile__) */
