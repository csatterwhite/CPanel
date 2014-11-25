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
    dXf /= rmax;
    
    // Do the same for derivative observation points
    for (int i=0; i<nb; i++)
    {
        dXb.row(i) = Xb.row(i)-X0;
    }
    dXb /= rmax;
    
    // Build Weights Matrices
    
    Eigen::VectorXd Wf(nf);
    Eigen::VectorXd Wb(nb);
    for (int i=0; i<nf; i++)
    {
        Wf(i) = 1/dXf.row(i).norm();
    }
    for (int i=0; i<nb; i++)
    {
        Wb(i) = 1/dXb.row(i).norm();
    }
    Wf.asDiagonal();
    Wb.asDiagonal();
    
    // Generate list of derivatives included in Taylor Series
    Eigen::MatrixXi ms = derivSequence(order,3);
    ms = sortBySum(ms);
    Eigen::VectorXi sums = ms.rowwise().sum();
    
    unsigned long t = ms.rows();
    
    // Build TLS and HTLS matrices (A and B) for function observations and derivative enforcement, respectively.
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nf,t);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nb,t);
    Eigen::Matrix<int,1,3> ps;
    double fact;
    double asum;
    for (int k=0; k<t; k++)
    {
        fact = 1/(factorial(ms(k,0))*factorial(ms(k,1))*factorial(ms(k,2)));
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
    
    // Build HTLS matrix for normal derivative enforcement
    if (!constr)
    {
        // Solve for unconstrained derivative coefficients
        Eigen::MatrixXd M = A.transpose()*Wf*A+B.transpose()*Wb*B;
        
        F = M.householderQr().solve(A.transpose()*Wf);
        G = M.householderQr().solve(B.transpose()*Wb);
        H = Eigen::MatrixXd::Zero(t,1);
    }
    else
    {
        // Solve for constrained derivative coefficents
        Eigen::MatrixXd V0full = Eigen::MatrixXd::Zero(t,1);
        V0full.block(0,0,2,0) = V0;
        
        // QR factorization
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(V0full);
        Eigen::MatrixXd Rfull = qr.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXd Q = qr.householderQ();
        
        double R = Rfull(0); // Only need first entrey of vector Rfull
        
        // Force R=1;
        Q /= R;
        R /= R;
        
        Eigen::MatrixXd M(A.rows()+B.rows()-1,A.cols());
        M << A,B;
        M *= Q;
        
        Eigen::MatrixXd Acon,Bcon,Alsq,Blsq;
        Acon = M.block(0, 0, nf-1, 0);
        Bcon = M.block(nf, 0, M.rows()-1, 0);
        Alsq = M.block(0, 1, nf-1, M.cols()-1);
        Blsq = M.block(nf, 1, M.rows()-1, M.cols()-1);
        
        M = Alsq.transpose()*Wf*Alsq+Blsq.transpose()*Wb*Blsq;
        
        Eigen::MatrixXd Fprm,Gprm,Hprm;
        Fprm = M.householderQr().solve(Alsq.transpose()*Wf);
        Gprm = M.householderQr().solve(Blsq.transpose()*Wb);
        Hprm = -(Fprm*Acon+Gprm*Bcon);
        
        Eigen::VectorXd Znf = Eigen::VectorXd::Zero(nf);
        Eigen::VectorXd Znb = Eigen::VectorXd::Zero(nb);
        F.resize(2,nf);
        G.resize(2,nb);
        H.resize(2,1);
        
        F << Znf,Fprm;
        G << Znb,Gprm;
        H << 1,Hprm;
        
        F = Q*F;
        G = Q*G;
        H = 1/R*Q*H;
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
    
    F /= rmax;
    G /= rmax;
    H /= rmax;
}

Eigen::MatrixXi chtlsnd::derivSequence(int q, int N)
{
    Eigen::MatrixXi d;
    if (q == 0)
    {
        d = Eigen::MatrixXi::Zero(1,N);
    }
    else
    {
        Eigen::MatrixXi sqenext;
        Eigen::VectorXi temp(N);
        for (int i=0; i<q; i++)
        {
            if (N > 1)
            {
                sqenext = derivSequence(q-i,N-1);
                for (int j=0; j<sqenext.rows(); j++)
                {
                    for (int k=0; k<sqenext.cols(); k++)
                    {
                        temp(k) = sqenext(j,k);
                    }
                    temp(N-1) = i;
                    d.row(d.rows()) = temp;
                }
            }
            else
            {
                d(d.rows(),1) = i;
            }
        }
    }
    return d;
}

Eigen::MatrixXi chtlsnd::sortBySum(Eigen::MatrixXi m)
{
    unsigned long rows = m.rows();
    unsigned long cols = m.cols();
    Eigen::MatrixXi out(rows,cols);
    Eigen::MatrixXi temp(rows,cols);
    out.resize(1,rows);
    out.row(0) = m.row(0);
    for (int i=1; i<m.rows(); i++)
    {
        for (int j=0; j<out.rows(); j++)
        {
            if (m.row(i).sum() < out.row(j).sum())
            {
                out = insertRow(out,m.row(i),j);
            }
        }
    }
    return out;
}

Eigen::MatrixXi chtlsnd::insertRow(const Eigen::Matrix<int,Eigen::Dynamic,3> &m, const Eigen::Matrix<int,1,3> &insert, int row)
{
    Eigen::Matrix<int,Eigen::Dynamic,3> temp(m.rows()+1,m.cols());
    unsigned long rows = m.rows();
    unsigned long cols = m.cols();
    assert(row < rows);
    if (row == 0)
    {
        temp.row(row) = insert;
        temp.block(1,0,rows,cols-1) = m.block(0,0,rows-1,cols-1);
    }
    else if (row == rows)
    {
        temp = m;
        temp.row(row) = insert;
    }
    else
    {
        temp.row(row) = insert;
        temp.block(0,0,row-1,m.cols()) = m.block(0,0,row-1,m.cols());
        temp.block(row+1, 0, rows, cols-1) = m.block(row,0,rows-1,cols-1);
    }
    return temp;
}

int chtlsnd::factorial(int i)
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