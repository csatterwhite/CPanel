//
//  main.cpp
//  WriteOutputFiles
//
//  Created by Chris Satterwhite on 11/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "geometry.h"
#include "octreeFile.h"

int main(int argc, const char * argv[])
{
    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
    std::string inFile = "genericAC_noWake.tri";

    geometry temp(path+inFile);
    octreeFile file(path + "genericAC_octree.txt",temp.getOctree());
    
}
