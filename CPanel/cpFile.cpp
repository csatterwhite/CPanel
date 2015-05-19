//
//  cpFile.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/28/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "cpFile.h"

cpFile::cpFile(std::string filename) : file(filename)
{
    parsefile();
}

cpFile::cpFile(const char* filename)
{
    file = filename;
    parsefile();
}

void cpFile::parsefile()
{
    std::size_t pathEnd = file.find_last_of("/")+1;
    std::size_t nameEnd = file.find_last_of(".");
    bool relPath = false;
    if (file.substr(0,1) != "/")
    {
        relPath = true;
    }
    // Get extension
    ext = file.substr(nameEnd,file.size()-nameEnd);
    name = file.substr(pathEnd,nameEnd-pathEnd);
    
    if (!relPath)
    {
        path = file.substr(0,pathEnd);
    }
    else
    {
        std::string relativePath = file.substr(0,pathEnd);
        std::string currentPath = boost::filesystem::current_path().string();
        
        std::stringstream ss;
        ss << currentPath << "/" << relativePath;
        path = ss.str();
        std::stringstream ss2;
        ss2 << path << name << ext;
        file = ss2.str();
    }
    
//    // Set geometry path and geometry name
//    if (pathEnd != 0)
//    {
//        path = file.substr(0,pathEnd);
//    }
//    else
//    {
//        boost::filesystem::path p = boost::filesystem::current_path();
//        name = file.substr(0,nameEnd);
//        changePath(p.string());
//    }
}

void cpFile::changePath(std::string newPath)
{
    path = newPath;
    std::stringstream temp;
    temp << newPath << "/" << name << ext;
    file = temp.str();
}