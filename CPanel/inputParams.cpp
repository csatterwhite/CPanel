//
//  inputParams.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/16/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "inputParams.h"

bool inputParams::set()
{
    std::ifstream fid;
    fid.open(inputFile->file);
    if (fid)
    {
        // Read in parameters from input file
        std::string str1,str2;
        std::string geomF;
        int n;
        fid >> str1;
        while (str1 != "GeomFile")
        {
            fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            fid >> str1;
            if (fid.eof())
            {
                std::cout << "ERROR : Input file could not be read" << std::endl;
                return false;
            }
        }
        
        fid >> str2 >> geomF;
        geomFile = new cpFile(geomF);
        if (!checkGeomFile())
        {
            return false;
        }
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
        makeWorkingDir();
        return true;
    }
    else
    {
        std::cout << "ERROR : Input file could not be opened" << std::endl;
        return false;
    }
}

bool inputParams::checkGeomFile()
{
    if (geomFile->ext == ".tri")
    {
        normFlag = false;
    }
    else if (geomFile->ext == ".tricp")
    {
        normFlag = true;
    }
    else
    {
        std::cout << "ERROR : Unsupported File Type \nAccepted Filetypes : '.tri','.tricp'" << std::endl;
        return false;
    }
    return true;
}

void inputParams::print(std::ostream &stream)
{
    int nChars = 18;
    stream << std::setw(nChars) << "Geometry File " << "-> " << geomFile->file << std::endl;
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

void inputParams::makeWorkingDir()
{
    std::string inPath = inputFile->path;
    inPath = inPath.substr(0,inPath.size()-1); // Remove trailing slash;
    
    std::size_t folderStart = inPath.find_last_of("/")+1;
    std::string inFolder = inPath.substr(folderStart,inPath.size()-folderStart);
    
    if (inFolder != inputFile->name)
    {
        std::stringstream subdir;
        subdir << inputFile->path << inputFile->name;
        boost::filesystem::path p = subdir.str();
        if (!boost::filesystem::exists(p))
        {
            boost::filesystem::create_directories(p);
        }
        chdir(subdir.str().c_str());
        boost::filesystem::path oldPath = geomFile->file;
        geomFile->changePath(boost::filesystem::current_path().string());
        boost::filesystem::path newPath = geomFile->file;
        boost::filesystem::rename(oldPath,newPath);
        writeInputFile();
        boost::filesystem::remove(inputFile->file);
    }
    else if (boost::filesystem::current_path() != inPath)
    {
        chdir(inPath.c_str());
    }
}

void inputParams::writeInputFile()
{
    std::ofstream fid;
    std::stringstream newInFile;
    newInFile << inputFile->name << inputFile->ext;
    fid.open(newInFile.str());
    
    fid << "%% CPanel Input File %%" << std::endl;
    fid << "GeomFile =\t" << geomFile->file << std::endl;
    fid << "Sref =\t" << Sref << std::endl;
    fid << "bref =\t" << bref << std::endl;
    fid << "cref =\t" << cref << std::endl;
    fid << "X_cg =\t" << cg(0) << std::endl;
    fid << "Y_cg =\t" << cg(1) << std::endl;
    fid << "Z_cg =\t" << cg(2) << std::endl;
    fid << "Velocity (ft/s)" << std::endl;
    fid << velocities.size() << std::endl;
    for (int i=0; i<velocities.size(); i++)
    {
        fid << velocities(i) << std::endl;
    }
    fid << "Angle of Attack (degrees)" << std::endl;
    fid << alphas.size() << std::endl;
    for (int i=0; i<alphas.size(); i++)
    {
        fid << alphas(i) << std::endl;
    }
    fid << "Angle of Sideslip (degrees)" << std::endl;
    fid << betas.size() << std::endl;
    for (int i=0; i<betas.size(); i++)
    {
        fid << betas(i) << std::endl;
    }
    fid << "Mach Number" << std::endl;
    fid << machs.size() << std::endl;
    for (int i=0; i<machs.size(); i++)
    {
        fid << machs(i) << std::endl;
    }
    fid.close();
}