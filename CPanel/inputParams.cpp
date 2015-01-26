//
//  inputParams.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/16/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "inputParams.h"

void inputParams::set(std::ifstream &fid)
{
    // Read in parameters from input file
    std::string str1,str2;
    int n;
    fid >> str1 >> str2 >> geomFile;
    checkGeomFile();
    fid >> str1 >> str2 >> Sref;
    fid >> str1 >> str2 >> bref;
    fid >> str1 >> str2 >> cref;
    for (int i=0; i<3; i++)
    {
        fid >> str1 >> str2 >> cg(i);
    }
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid >> n;
    velocities.resize(n);
    for (int i=0; i<n; i++)
    {
        fid >> velocities(i);
    }
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid >> n;
    alphas.resize(n);
    for (int i=0; i<n; i++)
    {
        fid >> alphas(i);
    }
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid >> n;
    betas.resize(n);
    for (int i=0; i<n; i++)
    {
        fid >> betas(i);
    }
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    fid >> n;
    machs.resize(n);
    for (int i=0; i<n; i++)
    {
        fid >> machs(i);
    }
}

void inputParams::checkGeomFile()
{
    
    std::size_t pathEnd = geomFile.find_last_of("/")+1;
    std::size_t nameEnd = geomFile.find_last_of(".");
    // Check file type
    std::string ext = geomFile.substr(nameEnd,geomFile.size()-nameEnd);
    if (ext == ".tri")
    {
        normFlag = false;
    }
    else if (ext == ".tricp")
    {
        normFlag = true;
    }
    else
    {
        std::cout << "ERROR : Unsupported File Type \nAccepted Filetypes : '.tri','.tricp'" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Set geometry path and geometry name
    if (pathEnd != 0)
    {
        geomPath = geomFile.substr(0,pathEnd);
        geomName = geomFile.substr(pathEnd,nameEnd-pathEnd);
    }
    else
    {
        boost::filesystem::path p = boost::filesystem::current_path();
        geomPath = p.string()+"/";
        geomName = geomFile.substr(0,nameEnd);
    }
    
    pathEnd = inputFile.find_last_of("/")+1;
    nameEnd = inputFile.find_last_of(".");
    
    // Set input path and input name
    if (pathEnd != 0)
    {
        inputPath = geomFile.substr(0,pathEnd);
        inputName = inputFile.substr(pathEnd,nameEnd-pathEnd);
    }
    else
    {
        boost::filesystem::path p = boost::filesystem::current_path();
        inputPath = p.string()+"/";
        inputName = inputFile.substr(0,nameEnd);
    }
}

void inputParams::print(std::ostream &stream)
{
    int nChars = 18;
    stream << std::setw(nChars) << "Geometry File " << "-> " << geomFile << std::endl;
    stream << std::setw(nChars) << "Reference Area " << "-> " << Sref << " ft^2" << std::endl;
    stream << std::setw(nChars) << "Reference Span " << "-> " << bref << " ft" << std::endl;
    stream << std::setw(nChars) << "Reference Chord " << "-> " << cref << " ft" << std::endl;
    stream << std::setw(nChars) << "Center of Gravity " << "-> [" << cg(0) << " " << cg(1) << " " << cg(2) << "] ft" << std::endl;
    
    stream << std::setw(nChars) << "Velocity " << "-> ";
    printVec(velocities,stream);
    stream << "ft/s" << std::endl;
    
    stream << std::setw(nChars) << "Alpha " << "-> ";
    printVec(alphas,stream);
    stream << "degrees" << std::endl;
    
    stream << std::setw(nChars) << "Beta " << "-> ";
    printVec(betas,stream);
    stream << "degrees" << std::endl;
    
    stream << std::setw(nChars) << "Mach # " << "-> ";
    printVec(machs,stream);
    stream << std::endl;
}

void inputParams::printVec(Eigen::VectorXd &vec,std::ostream &stream)
{
    if (vec.size() > 1)
    {
        for (int i=0; i<vec.size(); i++)
        {
            if (i==0)
            {
                stream << "[" << vec(i) << " ";
            }
            else if (i==vec.size()-1)
            {
                stream << vec(i) << "] ";
            }
            else
            {
                stream << vec(i) << " ";
            }
        }
    }
    else
    {
        stream << vec(0) << " ";
    }
}