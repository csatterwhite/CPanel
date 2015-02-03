//
//  cpFile.h
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/28/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__cpFile__
#define __CPanel___Unstructured_Panel_Code__cpFile__

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

struct cpFile
{
    std::string file;
    std::string path;
    std::string name;
    std::string ext;
    
    cpFile(std::string filename);
    cpFile(const char* filename);
    
    void changePath(std::string newPath);
    
private:
    void parsefile();
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__cpFile__) */
