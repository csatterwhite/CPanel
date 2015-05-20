//
//  CPanelMgr.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/16/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "CPanelMgr.h"

void caseMgr::setCases()
{
    cpCase* c;
    for (int v=0; v<p->velocities.rows(); v++)
    {
        for (int a=0; a<p->alphas.rows(); a++)
        {
            for (int b=0; b<p->betas.rows(); b++)
            {
                for (int m=0; m<p->machs.rows(); m++)
                {
                            c = new cpCase(geom,p->velocities(v),p->alphas(a),p->betas(b),p->machs(m),p);
                            cases.push_back(c);
                }
            }
        }
    }
}

void caseMgr::runCases()
{
    std::cout << "\nRunning " << cases.size() << " Cases... (\u2713 - Complete, X - Not Requested)\n" << std::endl;
    std::cout << std::setw(10) << std::left << "Case #" << std::setw(15) << std::left << "Solve System" << std::setw(15) << std::left << "Surface Data" << std::setw(16) << std::left << "Trefftz Plane" <<  std::setw(14) << std::left << "Streamlines" << std::setw(22) << std::left << "Stability Derivatives" << std::endl;
    for (int i=0; i<cases.size(); i++)
    {
        std::string out;
        std::stringstream outstream;
        outstream << i+1 << "/" << cases.size();
        out = outstream.str();
        std::cout << std::setw(10) << std::left << out << std::flush;
        cases[i]->run(true,p->surfStreamFlag,p->stabDerivFlag);
    }
    std::cout << "Complete" << std::endl;
}

void caseMgr::writeSummary()
{
    std::ofstream out;
    std::string outFile = p->inputFile->name+".CPout";
    out.open(outFile);
    if (out)
    {
        out << "-----------INPUTS-----------\n\n";
        p->print(out);
        out << "\n\n";
        out << "-----------RESULTS-----------\n\n";
        out << std::left; // Left justify results
        out << std::fixed;
        out << std::setprecision(3);
        outSpacing.resize(9);
        outSpacing << 11,7,8,8,8,8,12,12,12;
        for (int i=0; i<cases.size(); i++)
        {
            writeCase(i+1, cases[i], out);
        }
    }
    std::cout << "Case data written to " << outFile << std::endl;
    
}

void caseMgr::writeCase(int caseNumber, cpCase* c, std::ofstream &outStream)
{
    outStream << std::endl;
    outStream << "---Case #" << caseNumber << "---" << std::endl;
    outStream << "\n\t--Flow Conditions--\n" << std::endl;
    outStream << "\t\t" << std::setw(outSpacing(0)) << "Velocity" << std::setw(outSpacing(1)) << "Mach" << std::setw(outSpacing(2)) << "Alpha" << std::setw(outSpacing(3)) << "Beta" << "\n";
    outStream << "\t\t" << std::setw(outSpacing(0)) << c->getV() << std::setw(outSpacing(1)) << c->getMach() << std::setw(outSpacing(2)) << c->getAlpha() << std::setw(outSpacing(3)) << c->getBeta() << "\n";
    outStream << "\n\t--Force Coefficients--\n" << std::endl;
    outStream << "\t\t-Trefftz Plane-" << std::endl;
    outStream << "\t\t\tCL = " << c->getCL() << "\tCDi = " << c->getCD() << std::endl;
    outStream << "\n\t\t-Surface Integrated Force Coefficients-" << std::endl;
    outStream << "\t\t\t" << std::setw(15) << "Body Axis" << std::setw(8) << "CN" << std::setw(8) << "CA" << std::setw(8) << "CY" << std::endl;
    outStream << "\t\t\t" << std::setw(15) << " " << std::setw(8) << c->getBodyForces()(2) << std::setw(8) << c->getBodyForces()(0) << std::setw(8) << c->getBodyForces()(1) << std::endl;
    outStream << "\t\t\t" << std::setw(15) << "Wind Axis" << std::setw(8) << "CL" << std::setw(8) << "CD" << std::setw(8) << "CY" << std::endl;
    outStream << "\t\t\t" << std::setw(15) << " " << std::setw(8) << c->getWindForces()(2) << std::setw(8) << c->getWindForces()(0) << std::setw(8) << c->getWindForces()(1) << std::endl;
    
    outStream << "\n\t--Moment Coefficients--\n" << std::endl;
    outStream << "\t\t" << std::setw(12) << "Cm (pitch)" << std::setw(12) << "Cl (roll)" << std::setw(12) << "Cn (yaw)" << std::endl;
    outStream << "\t\t" << std::setw(12) << c->getMoment()(1) << std::setw(12) << c->getMoment()(0) << std::setw(12) << c->getMoment()(2) << std::endl;
    
    outStream << "\n\t--Stability Derivatives--\n" << std::endl;
    outStream << "\t\t" << std::setw(12) << "CL_alpha" << std::setw(12) << "CY_alpha" << std::setw(12) << "Cm_alpha" << std::setw(12) << "Cl_alpha" << std::setw(12) << "Cn_alpha" << std::endl;
    outStream << "\t\t" << std::setw(12) << c->get_dF_dAlpha()(2) << std::setw(12) << c->get_dF_dAlpha()(1) << std::setw(12) << c->get_dM_dAlpha()(1) << std::setw(12) << c->get_dM_dAlpha()(0) << std::setw(12) << c->get_dM_dAlpha()(2) << std::endl;
    
    outStream << "\t\t" << std::setw(12) << "CL_beta" << std::setw(12) << "CY_beta" << std::setw(12) << "Cm_beta" << std::setw(12) << "Cl_beta" << std::setw(12) << "Cn_beta" << std::endl;
    outStream << "\t\t" << std::setw(12) << c->get_dF_dBeta()(2) << std::setw(12) << c->get_dF_dBeta()(1) << std::setw(12) << c->get_dM_dBeta()(1) << std::setw(12) << c->get_dM_dBeta()(0) << std::setw(12) << c->get_dM_dBeta()(2) << std::endl;
    outStream << "\n\n" << std::endl;
    
}
