//
//  chtlsnd.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 11/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "chtlsnd.h"

chtlsnd::chtlsnd(const Eigen::Matrix<double,1,Eigen::Dynamic> &X0, const Eigen::MatrixXd &Xf, int order, const Eigen::MatrixXd &Xb, const Eigen::MatrixXd &Vb, Eigen::VectorXd V0)
{
    unsigned long nf = Xf.rows();
    unsigned long nb = Xb.rows();
    
    
    bool constr = true;
    if (V0.isZero(0))
    {
        constr = false;
    }
    int N = (int)X0.size(); // Problem Dimensionality
    
    Eigen::MatrixXd dXf(nf,N);
    Eigen::MatrixXd dXb(nb,N);
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
    
    Eigen::VectorXd Wf = Eigen::VectorXd::Zero(nf);
    Eigen::VectorXd Wb = Eigen::VectorXd::Zero(nb);
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
    
    // Generate list of derivatives included in Taylor Series
    Eigen::MatrixXi ms = derivSequence(order,N);
    ms = sortBySum(ms);
    ms = ms.block(1, 0, ms.rows()-1, ms.cols());
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
        fact = 1;
        for (int i=0; i<N; i++)
        {
            fact /= factorial(ms(k,i));
        }
        for (int i=0; i<nf; i++)
        {
            double prod = 1;
            for (int ii=0; ii<N; ii++)
            {
                prod *= pow(dXf(i,ii),ms(k,ii));
            }
            A(i,k) = fact*prod;
//            A(i,k) = fact*(pow(dXf(i,0),ms(k,0))*pow(dXf(i,1),ms(k,1))*pow(dXf(i,2),ms(k,2)));
        }
        for (int i=0; i<nb; i++)
        {
            asum = 0;
            for (int a=0; a<N; a++)
            {
                if (ms(k,a) != 0) // Skip ^0 cases
                {
                    ps = ms.row(k);
                    ps(a) -= 1; // Differentiat by loweing k'th power
                    double prod = 1;
                    for (int ii=0; ii<N; ii++)
                    {
                        prod *= pow(dXb(i,ii),ps(ii));
                    }
                    asum += Vb(i,a)*ms(k,a)*prod;
//                    asum += Vb(i,a)*ms(k,a)*(pow(dXb(i,0),ps(0))*pow(dXb(i,1),ps(1))*pow(dXb(i,2),ps(2)));
                }
            }
            B(i,k) = fact*asum;
        }
    }

    // Build HTLS matrix for normal derivative enforcement
    if (!constr)
    {
        // Solve for unconstrained derivative coefficients
        Eigen::MatrixXd M(A.rows()-1,A.cols());
        M = A.transpose()*WfDiag*A+B.transpose()*WbDiag*B;
        
        F = M.householderQr().solve(A.transpose()*WfDiag);
        G = M.householderQr().solve(B.transpose()*WbDiag);
        H = Eigen::MatrixXd::Zero(t,1);
    }
    else
    {
        // Solve for constrained derivative coefficents
        Eigen::MatrixXd V0full = Eigen::MatrixXd::Zero(t,1);
        V0full.block(0,0,N,1) = V0;

        // QR factorization
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(V0full);
        Eigen::MatrixXd Rfull = qr.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXd Q = qr.householderQ();

        double R = Rfull(0); // Only need first entry of vector Rfull
      
        // Force R=1;
        Q /= R;
        R /= R;

        Eigen::MatrixXd M(A.rows()+B.rows(),A.cols());
        M << A,B;
        M *= Q;

        Eigen::MatrixXd Acon,Bcon,Alsq,Blsq;
        Acon = M.block(0, 0, nf, 1);
        Bcon = M.block(nf, 0, M.rows()-nf, 1);
        Alsq = M.block(0, 1, nf, M.cols()-1);
        Blsq = M.block(nf, 1, M.rows()-nf, M.cols()-1);

        M = Alsq.transpose()*WfDiag*Alsq+Blsq.transpose()*WbDiag*Blsq;
        
        Eigen::MatrixXd Fprm,Gprm,Hprm;
        Fprm = M.householderQr().solve(Alsq.transpose()*WfDiag);
        Gprm = M.householderQr().solve(Blsq.transpose()*WbDiag);
        Hprm = -(Fprm*Acon + Gprm*Bcon);

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

}

Eigen::MatrixXi chtlsnd::derivSequence(int q, int N)
{
    // Builds the list of mixed partial derivatives to highest order q for a function of N dimensions. The algorithm used is recursive.
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
            out = insertRow(out,m.row(i),(int)out.rows());
        }
    }
    return out;
    
}

Eigen::MatrixXi chtlsnd::insertRow(const Eigen::MatrixXi &m, const Eigen::MatrixXi &insert, int row)
{
    Eigen::MatrixXi temp(m.rows()+1,m.cols());
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