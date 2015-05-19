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
    // Set Default Parameters;
    Sref = 1;
    bref = 1;
    cref = 1;
    cg << 0,0,0;
    velocities.resize(1);
    alphas.resize(1);
    betas.resize(1);
    machs.resize(1);
    velocities(0) = 1;
    alphas(0) = 0;
    betas(0) = 0;
    machs(0) = 0.1;
    
    // Read parameters from input file
    std::ifstream fid;
    fid.open(inputFile->file);
    if (fid)
    {
        std::string s1,s2,s3;
        int n;
        std::getline(fid,s1);
        while (!fid.eof())
        {
            fid >> s1;
            if (s1.compare("%") == 0)
            {
                // Skip Comment Lines
                fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }
            else
            {
                if (s1.compare("GeomFile") == 0)
                {
                    fid >> s2 >> s3;
                    geomFile = new cpFile(s3);
                    checkGeomFile();
                }
                else if (s1.compare("S_ref") == 0)
                {
                    fid >> s2 >> Sref;
                }
                else if (s1.compare("b_ref") == 0)
                {
                    fid >> s2 >> bref;
                }
                else if (s1.compare("c_ref") == 0)
                {
                    fid >> s2 >> cref;
                }
                else if (s1.compare("X_cg") == 0)
                {
                    fid >> s2 >> cg(0);
                }
                else if (s1.compare("Y_cg") == 0)
                {
                    fid >> s2 >> cg(1);
                }
                else if (s1.compare("Z_cg") == 0)
                {
                    fid >> s2 >> cg(2);
                }
                else if (s1.compare("Velocity") == 0)
                {
                    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
                    fid >> n;
                    velocities.resize(n);
                    for (int i=0; i<n; i++)
                    {
                        fid >> velocities(i);
                    }
                }
                else if (s1.compare("Angle_of_Attack") == 0)
                {
                    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
                    fid >> n;
                    alphas.resize(n);
                    for (int i=0; i<n; i++)
                    {
                        fid >> alphas(i);
                    }
                }
                else if (s1.compare("Angle_of_Sideslip") == 0)
                {
                    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
                    fid >> n;
                    betas.resize(n);
                    for (int i=0; i<n; i++)
                    {
                        fid >> betas(i);
                    }
                }
                else if (s1.compare("Mach_Number") == 0)
                {
                    fid.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
                    fid >> n;
                    machs.resize(n);
                    for (int i=0; i<n; i++)
                    {
                        fid >> machs(i);
                    }
                }
                else if (s1.compare("Surface_Streamlines") == 0)
                {
                    fid >> surfStreamFlag;
                }
                else if (s1.compare("Stability_Derivatives") == 0)
                {
                    fid >> stabDerivFlag;
                }
                else if (s1.compare("Write_Influence_Coefficients") == 0)
                {
                    fid >> writeCoeffFlag;
                }
            }
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
    int nChars = 28;
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
    
    stream << std::setw(nChars) << "Surface Streamlines " << "-> ";
    if (surfStreamFlag)
        stream << "ON" << std::endl;
    else
        stream << "OFF" << std::endl;
    
    stream << std::setw(nChars) << "Stability Derivatives " << "-> ";
    if (stabDerivFlag)
        stream << "ON" << std::endl;
    else
        stream << "OFF" << std::endl;
    
    stream << std::setw(nChars) << "Write Inf Coeff to File " << "-> ";
    if (writeCoeffFlag)
        stream << "ON" << std::endl;
    else
        stream << "OFF" << std::endl;
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
    
    fid << "%% CPanel Input File %%\n" << std::endl;
    fid << "% Reference Geometry %" << std::endl;
    fid << "GeomFile =\t" << geomFile->file << std::endl;
    fid << "S_ref =\t" << Sref << std::endl;
    fid << "b_ref =\t" << bref << std::endl;
    fid << "c_ref =\t" << cref << std::endl;
    fid << "X_cg =\t" << cg(0) << std::endl;
    fid << "Y_cg =\t" << cg(1) << std::endl;
    fid << "Z_cg =\t" << cg(2) << std::endl;
    fid << std::endl;
    fid << "% Cases %" << std::endl;
    fid << "Velocity (ft/s)" << std::endl;
    fid << velocities.size() << std::endl;
    for (int i=0; i<velocities.size(); i++)
    {
        fid << velocities(i) << std::endl;
    }
    fid << "Angle_of_Attack (degrees)" << std::endl;
    fid << alphas.size() << std::endl;
    for (int i=0; i<alphas.size(); i++)
    {
        fid << alphas(i) << std::endl;
    }
    fid << "Angle_of_Sideslip (degrees)" << std::endl;
    fid << betas.size() << std::endl;
    for (int i=0; i<betas.size(); i++)
    {
        fid << betas(i) << std::endl;
    }
    fid << "Mach_Number" << std::endl;
    fid << machs.size() << std::endl;
    for (int i=0; i<machs.size(); i++)
    {
        fid << machs(i) << std::endl;
    }
    fid << std::endl;
    fid << "% Solver Options (0 = OFF, 1 = ON) %" << std::endl;
    fid << "Surface_Streamlines" << std::endl;
    fid << surfStreamFlag << std::endl;
    fid << "Stability_Derivatives" << std::endl;
    fid << stabDerivFlag << std::endl;
    fid << "Write_Influence_Coefficients" << std::endl;
    fid << writeCoeffFlag << std::endl;
    
    fid.close();
}