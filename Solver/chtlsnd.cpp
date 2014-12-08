//
//  chtlsnd.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 11/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "chtlsnd.h"

chtlsnd::chtlsnd(const Eigen::Matrix<double,1,3> &X0, const Eigen::Matrix<double,Eigen::Dynamic,3> &Xf, int order, const Eigen::Matrix<double,Eigen::Dynamic,3> &Xb, const Eigen::Matrix<double,Eigen::Dynamic,3> &Vb, Eigen::Vector3d V0)
{
    unsigned long nf = Xf.rows();
    unsigned long nb = Xb.rows();
    
    bool constr = true;
    if (V0.isZero())
    {
        constr = false;
    }
    
    Eigen::MatrixXd dXf(nf,3);
    Eigen::MatrixXd dXb(nb,3);
    double rmax = 0;
    
//    ////////////////////////
//    for (int i=0; i<Xf.rows(); i++)
//    {
//        for (int j=0; j<Xf.cols(); j++)
//        {
//            std::cout << Xf(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
//    ////////////////////////
//    for (int i=0; i<Xb.rows(); i++)
//    {
//        for (int j=0; j<Xb.cols(); j++)
//        {
//            std::cout << Xb(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
    
    // Calculate and Normalize relative locations of function observation points
    for (int i=0; i<nf; i++)
    {
        dXf.row(i) = Xf.row(i)-X0;
        double r = dXf.row(i).norm();
        if (r > rmax)
        {
            rmax = r;
        }
    }
//    ////////////////////////
//    for (int i=0; i<dXf.rows(); i++)
//    {
//        for (int j=0; j<dXf.cols(); j++)
//        {
//            std::cout << dXf(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
//    ////////////////////////
//    std::cout << rmax << std::endl;
//    ////////////////////////
    
    dXf /= rmax;
    
    // Do the same for derivative observation points
    for (int i=0; i<nb; i++)
    {
        dXb.row(i) = Xb.row(i)-X0;
    }
    
//    ////////////////////////
//    for (int i=0; i<dXb.rows(); i++)
//    {
//        for (int j=0; j<dXb.cols(); j++)
//        {
//            std::cout << dXb(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
    
    dXb /= rmax;
    
//    ////////////////////////
//    for (int i=0; i<dXf.rows(); i++)
//    {
//        for (int j=0; j<dXf.cols(); j++)
//        {
//            std::cout << dXf(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
//    ////////////////////////
//    for (int i=0; i<dXb.rows(); i++)
//    {
//        for (int j=0; j<dXb.cols(); j++)
//        {
//            std::cout << dXb(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
    
    // Build Weights Matrices
    
    Eigen::VectorXd Wf(nf), Wb(nb);
    Eigen::MatrixXd WfDiag(nf,nf), WbDiag(nb,nb);
    for (int i=0; i<nf; i++)
    {
        Wf(i) = 1/dXf.row(i).norm();
    }
    for (int i=0; i<nb; i++)
    {
        Wb(i) = 1/dXb.row(i).norm();
    }
    WfDiag = Wf.asDiagonal();
    WbDiag = Wb.asDiagonal();
    
//    ////////////////////////
//    for (int i=0; i<WfDiag.rows(); i++)
//    {
//        for (int j=0; j<WfDiag.cols(); j++)
//        {
//            std::cout << WfDiag(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
//    ////////////////////////
//    for (int i=0; i<WbDiag.rows(); i++)
//    {
//        for (int j=0; j<WbDiag.cols(); j++)
//        {
//            std::cout << WbDiag(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
    
    // Generate list of derivatives included in Taylor Series
    Eigen::MatrixXi ms = derivSequence(order,3);
    ms = sortBySum(ms);
    ms = ms.block(1, 0, ms.rows()-1, ms.cols());
    Eigen::VectorXi sums = ms.rowwise().sum();
    
//    ////////////////////////
//    for (int i=0; i<ms.rows(); i++)
//    {
//        for (int j=0; j<ms.cols(); j++)
//        {
//            std::cout << ms(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
    
    unsigned long t = ms.rows();
    
//    ///////////////////////
//    std::cout << t << std::endl;
//    ///////////////////////

    
    // Build TLS and HTLS matrices (A and B) for function observations and derivative enforcement, respectively.
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nf,t);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nb,t);
    Eigen::Matrix<int,1,3> ps;
    double fact;
    double asum;
    for (int k=0; k<t; k++)
    {
        fact = 1;
        for (int i=0; i<3; i++)
        {
            fact /= factorial(ms(k,i));
        }
        for (int i=0; i<nf; i++)
        {
            A(i,k) = fact*(pow(dXf(i,0),ms(k,0))*pow(dXf(i,1),ms(k,1))*pow(dXf(i,2),ms(k,2)));
        }
        for (int i=0; i<nb; i++)
        {
            asum = 0;
            for (int a=0; a<3; a++)
            {
                if (ms(k,a) != 0) // Skip ^0 cases
                {
                    ps = ms.row(k);
                    ps(a) -= 1; // Differentiat by loweing k'th power
                    asum += Vb(i,a)*ms(k,a)*(pow(dXb(i,0),ps(0))*pow(dXb(i,1),ps(1))*pow(dXb(i,2),ps(2)));
                }
            }
            B(i,k) = fact*asum;
        }
    }
    
//    ////////////////////////
//    for (int i=0; i<A.rows(); i++)
//    {
//        for (int j=0; j<A.cols(); j++)
//        {
//            std::cout << A(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
//    ////////////////////////
//    for (int i=0; i<B.rows(); i++)
//    {
//        for (int j=0; j<B.cols(); j++)
//        {
//            std::cout << B(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
    
    // Build HTLS matrix for normal derivative enforcement
    if (!constr)
    {
        // Solve for unconstrained derivative coefficients
        Eigen::MatrixXd M = A.transpose()*Wf*A+B.transpose()*Wb*B;
        
        F = M.householderQr().solve(A.transpose()*WfDiag);
        G = M.householderQr().solve(B.transpose()*WbDiag);
        H = Eigen::MatrixXd::Zero(t,1);
    }
    else
    {
        // Solve for constrained derivative coefficents
        Eigen::MatrixXd V0full = Eigen::MatrixXd::Zero(t,1);
        V0full.block(0,0,3,1) = V0;
        
//        ////////////////////////
//        for (int i=0; i<V0full.rows(); i++)
//        {
//            for (int j=0; j<V0full.cols(); j++)
//            {
//                std::cout << V0full(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
        
        // QR factorization
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(V0full);
        Eigen::MatrixXd Rfull = qr.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXd Q = qr.householderQ();
        
//        ////////////////////////
//        for (int i=0; i<Q.rows(); i++)
//        {
//            for (int j=0; j<Q.cols(); j++)
//            {
//                std::cout << Q(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
        
        double R = Rfull(0); // Only need first entry of vector Rfull
        
//        ////////////////////////
//        std::cout << R << std::endl;
//        ////////////////////////
        
        // Force R=1;
        Q /= R;
        R /= R;
        
//        ////////////////////////
//        for (int i=0; i<Q.rows(); i++)
//        {
//            for (int j=0; j<Q.cols(); j++)
//            {
//                std::cout << Q(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        std::cout << R << std::endl;
//        ////////////////////////
        
        Eigen::MatrixXd M(A.rows()+B.rows(),A.cols());
        M << A,B;
        M *= Q;
        
//        ////////////////////////
//        for (int i=0; i<M.rows(); i++)
//        {
//            for (int j=0; j<M.cols(); j++)
//            {
//                std::cout << M(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
        
        Eigen::MatrixXd Acon,Bcon,Alsq,Blsq;
        Acon = M.block(0, 0, nf, 1);
        Bcon = M.block(nf, 0, M.rows()-nf, 1);
        Alsq = M.block(0, 1, nf, M.cols()-1);
        Blsq = M.block(nf, 1, M.rows()-nf, M.cols()-1);
        
//        ////////////////////////
//        for (int i=0; i<Acon.rows(); i++)
//        {
//            for (int j=0; j<Acon.cols(); j++)
//            {
//                std::cout << Acon(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        for (int i=0; i<Bcon.rows(); i++)
//        {
//            for (int j=0; j<Bcon.cols(); j++)
//            {
//                std::cout << Bcon(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        for (int i=0; i<Alsq.rows(); i++)
//        {
//            for (int j=0; j<Alsq.cols(); j++)
//            {
//                std::cout << Alsq(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        for (int i=0; i<Blsq.rows(); i++)
//        {
//            for (int j=0; j<Blsq.cols(); j++)
//            {
//                std::cout << Blsq(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
        
        M = Alsq.transpose()*WfDiag*Alsq+Blsq.transpose()*WbDiag*Blsq;
//        ////////////////////////
//        for (int i=0; i<M.rows(); i++)
//        {
//            for (int j=0; j<M.cols(); j++)
//            {
//                std::cout << M(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
        
        Eigen::MatrixXd Fprm,Gprm,Hprm;
        Fprm = M.householderQr().solve(Alsq.transpose()*WfDiag);
        Gprm = M.householderQr().solve(Blsq.transpose()*WbDiag);
        Hprm = -(Fprm*Acon+Gprm*Bcon);
        
//        ////////////////////////
//        for (int i=0; i<Fprm.rows(); i++)
//        {
//            for (int j=0; j<Fprm.cols(); j++)
//            {
//                std::cout << Fprm(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        for (int i=0; i<Gprm.rows(); i++)
//        {
//            for (int j=0; j<Gprm.cols(); j++)
//            {
//                std::cout << Gprm(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        for (int i=0; i<Hprm.rows(); i++)
//        {
//            for (int j=0; j<Hprm.cols(); j++)
//            {
//                std::cout << Hprm(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
        
        Eigen::MatrixXd Znf = Eigen::MatrixXd::Zero(1,nf);
        Eigen::MatrixXd Znb = Eigen::MatrixXd::Zero(1,nb);
        F.resize(Fprm.rows()+1,nf);
        G.resize(Gprm.rows()+1,nb);
        H.resize(Hprm.rows()+1,1);
        
        F << Znf,Fprm;
        G << Znb,Gprm;
        H << 1,Hprm;
        
        F = Q*F;
        G = Q*G;
        H = 1/R*Q*H;
        
//        ////////////////////////
//        for (int i=0; i<F.rows(); i++)
//        {
//            for (int j=0; j<F.cols(); j++)
//            {
//                std::cout << F(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        for (int i=0; i<G.rows(); i++)
//        {
//            for (int j=0; j<G.cols(); j++)
//            {
//                std::cout << G(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
//        ////////////////////////
//        for (int i=0; i<H.rows(); i++)
//        {
//            for (int j=0; j<H.cols(); j++)
//            {
//                std::cout << H(i,j) << "\t";
//            }
//            std::cout << std::endl;
//        }
//        ////////////////////////
    }
    // Redimensionalize derivative coefficients;
    for (int i=0; i<F.rows(); i++)
    {
        for (int j=0; j<F.cols(); j++)
        {
            F(i,j) /= pow(rmax,sums(i));
        }
    }
    for (int i=0; i<G.rows(); i++)
    {
        for (int j=0; j<G.cols(); j++)
        {
            G(i,j) /= pow(rmax,sums(i)-1);
        }
    }
    for (int i=0; i<H.rows(); i++)
    {
        H(i) = H(i)/pow(rmax,sums(i)-1);
    }
    
//    ////////////////////////
//    for (int i=0; i<F.rows(); i++)
//    {
//        for (int j=0; j<F.cols(); j++)
//        {
//            std::cout << F(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
//    ////////////////////////
//    for (int i=0; i<G.rows(); i++)
//    {
//        for (int j=0; j<G.cols(); j++)
//        {
//            std::cout << G(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
//    ////////////////////////
//    for (int i=0; i<H.rows(); i++)
//    {
//        for (int j=0; j<H.cols(); j++)
//        {
//            std::cout << H(i,j) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    ////////////////////////
}

Eigen::MatrixXi chtlsnd::derivSequence(int q, int N)
{
    Eigen::MatrixXi d(0,0);
    if (q == 0)
    {
        d = Eigen::MatrixXi::Zero(1,N);
    }
    else
    {
        Eigen::MatrixXi sqenext;
        Eigen::VectorXi temp;
        int count = 0;
        for (int i=0; i<=q; i++)
        {
            if (N > 1)
            {
                sqenext = derivSequence(q-i,N-1);
                for (int j=0; j<sqenext.rows(); j++)
                {
                    temp = sqenext.row(j);
                    temp.conservativeResize(N);
                    temp(N-1) = i;
                    d.conservativeResize(d.rows()+1,N);
                    d.row(d.rows()-1) = temp;
                    count++;
                }
            }
            else
            {
                d.conservativeResize(d.rows()+1,N);
                d(i,0) = i;
            }
        }
    }
    return d;
}

Eigen::MatrixXi chtlsnd::sortBySum(Eigen::MatrixXi m)
{
    unsigned long cols = m.cols();
    Eigen::MatrixXi out(1,cols);
    out.row(0) = m.row(0);
    bool flag;
    for (int i=1; i<m.rows(); i++)
    {
        flag = false;
        for (int j=0; j<out.rows(); j++)
        {
            if (m.row(i).sum() < out.row(j).sum())
            {
                out = insertRow(out,m.row(i),j);
                flag = true;
                break;
            }
        }
        if (!flag)
        {
            out = insertRow(out,m.row(i),out.rows());
        }
    }
    return out;
    
}

Eigen::MatrixXi chtlsnd::insertRow(const Eigen::Matrix<int,Eigen::Dynamic,3> &m, const Eigen::Matrix<int,1,3> &insert, int row)
{
    Eigen::Matrix<int,Eigen::Dynamic,3> temp(m.rows()+1,m.cols());
    unsigned long rows = m.rows();
    unsigned long cols = m.cols();
    assert(row <= rows);
    if (row == 0)
    {
        temp.row(row) = insert;
        temp.block(1,0,rows,cols) = m.block(0,0,rows,cols);
    }
    else if (row == rows)
    {
        temp = m;
        temp.conservativeResize(m.rows()+1,m.cols());
        temp.row(row) = insert;
    }
    else
    {
        temp.block(0,0,row,cols) = m.block(0,0,row,cols);
        temp.block(row+1, 0, rows-(row), cols) = m.block(row,0,rows-(row),cols);
        temp.row(row) = insert;
    }
    return temp;
}

int chtlsnd::factorial(int i)
{
    if (i==0)
    {
        return 1;
    }
    else
    {
        int out = i;
        int mult = i-1;
        while (mult != 0)
        {
            out *= mult;
            mult--;
        }
        return out;
    }
}