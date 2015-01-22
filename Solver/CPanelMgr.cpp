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
    for (int v=0; v<p.velocities.rows(); v++)
    {
        for (int a=0; a<p.alphas.rows(); a++)
        {
            for (int b=0; b<p.betas.rows(); b++)
            {
                for (int m=0; m<p.machs.rows(); m++)
                {
                            c = new cpCase(geom,p.velocities(v),p.alphas(a),p.betas(b),p.machs(m),p.inputPath,p.inputName);
                            cases.push_back(c);
                }
            }
        }
    }
}

void caseMgr::runCases()
{
    std::cout << "Running Cases..." << std::endl;
    std::cout << std::setw(10) << std::left << "Case #" << std::setw(17) << std::left << "Solve System" << std::setw(14) << std::left << "Flow Data" << std::setw(18) << std::left << "Trefftz Plane" << std::endl;
    for (int i=0; i<cases.size(); i++)
    {
        std::string out;
        std::stringstream outstream;
        outstream << i+1 << "/" << cases.size();
        out = outstream.str();
        std::cout << std::setw(10) << std::left << out << std::flush;
        cases[i]->run();
    }
    std::cout << "Complete" << std::endl;
}

void caseMgr::writeSummary()
{
    std::ofstream out;
    std::string outFile = p.inputPath+p.inputName+".CPout";
    out.open(outFile);
    if (out)
    {
        out << "-----------INPUTS-----------\n\n";
        p.print(out);
        out << "\n\n";
        out << "-----------RESULTS-----------\n\n";
        out << std::left; // Left justify results
        out << std::fixed;
        out << std::setprecision(3);
        Eigen::VectorXi s(9);
        s << 11,7,8,8,8,8,12,12,12;
        out << std::setw(s(0)) << "Velocity" << std::setw(s(1)) << "Mach" << std::setw(s(2)) << "Alpha" << std::setw(s(3)) << "Beta" << std::setw(s(4)) << "CL" << std::setw(s(5)) << "CD" << std::setw(s(6)) << "Cm (pitch)" << std::setw(s(7)) << "Cl (roll)" << std::setw(s(8)) << "Cn (yaw)" << "\n";
        for (int i=0; i<cases.size(); i++)
        {
            out << std::setw(s(0)) << cases[i]->getV() << std::setw(s(1)) << cases[i]->getMach() << std::setw(s(2)) << cases[i]->getAlpha() << std::setw(s(3)) << cases[i]->getBeta() << std::setw(s(4)) << cases[i]->getCL() << std::setw(s(5)) << cases[i]->getCD() << std::setw(s(6)) << cases[i]->getMoment()(1) << std::setw(s(7)) << cases[i]->getMoment()(0) << std::setw(s(8)) << cases[i]->getMoment()(2) << "\n";
        }
    }
    std::cout << "Case data written to " << outFile << std::endl;
    
}
