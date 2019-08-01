//
//  The SR algorithm for the MLE for the completely monotone distribution
//
//  Created by Piet Groeneboom on 6/8/2019.
//  Copyright (c) 2019 Piet Groeneboom. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

void MLE(int n, int m, int discr[], double freq[], double cumw[], double cs[],
         double zz[], double pp[], double F[]);
void concavemaj(int n, double cumw[], double cs[],  double y[]);
  

// [[Rcpp::export]]

List monotone_1(DataFrame input)
{
    int         i,m,n;
    int         *discr;
    double      *F,*pp,*freq;
    double      *cumw,*cs,*zz;
    
    
    DataFrame DF = Rcpp::DataFrame(input);
    IntegerVector xcoor = DF["V1"];
    
    n = (int)xcoor.size();
    
    discr = new int[n+1];
    
    pp =    new double[n+1];
    F =    new double[n+1];
    freq    = new double[n+1];
    
    cumw    = new double[n+1];
    cs      = new double[n+1];
    zz      = new double[n+1];
    
    
    for (i=0;i<n;i++)
        discr[i]=(double)xcoor[i];
    
    m=0;
    for (i=0;i<n;i++)
    {
        if (discr[i]>m)
            m=discr[i];
    }
    
    MLE(n,m,discr,freq,cumw,cs,zz,pp,F);
    
    IntegerVector out0 = IntegerVector(n);
    
    for (i=0;i<n;i++)
        out0(i)=discr[i];
    
    NumericMatrix out1 = NumericMatrix(m+1,2);
    
    for (i=0;i<=m;i++)
    {
        out1(i,0)=i;
        out1(i,1)=pp[i];
    }
    
    NumericMatrix out2 = NumericMatrix(m+1,2);
    
    for (i=0;i<=m;i++)
    {
        out2(i,0)=i;
        out2(i,1)=F[i];
    }
    
    double out3 = F[10];
    
    // make the list for the output, containing the two estimates and the log likelihood
    
    List monotone_out = List::create(Rcpp::Named("data")=out0,Rcpp::Named("MLE")=out1,Rcpp::Named("F")=out2,Rcpp::Named("value")=out3);

    // free memory
    
    delete[] discr;  delete[] pp; delete[] F; delete[] freq;
    delete[] cumw; delete[] cs; delete[] zz;

    return monotone_out;
    
}

void MLE(int n, int m, int discr[], double freq[], double cumw[], double cs[],  double zz[], double pp[], double F[])
{
    int i,j;
    
    for (i=0;i<=m;i++)
        freq[i]=0;
    
    for (i=0;i<=m;i++)
    {
        for (j=0;j<n;j++)
        {
            if (discr[j]==i)
                freq[i] += 1.0/n;
        }
    }
    
    for (i=1;i<=m;i++)
        freq[i] = freq[i-1] +freq[i];
    
    cumw[0]=cs[0]=0;
    for (i=1;i<=m+1;i++)
    {
        cumw[i]=i*1.0;
        cs[i]=freq[i-1];
    }
    
    concavemaj(m+1,cumw,cs,zz);
    
    for (i=1;i<=m;i++)
        pp[i]=i*(zz[i]-zz[i+1]);
    
    pp[m+1]=(m+1)*zz[m+1];
    
    zz[0]=0;
    
    for (i=1;i<=m+1;i++)
        zz[i]=zz[i-1]+pp[i];
    
    for (i=0;i<=m;i++)
        F[i]=zz[i+1];
}

void concavemaj(int n, double cumw[], double cs[],  double y[])
{
  int    i,j,m;
  
  y[1] = cs[1]/cumw[1];
  for (i=2;i<=n;i++)
  {
    if (cumw[i]-cumw[i-1]>0)
      y[i] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
    else
      y[i]=0;
    if (y[i-1]<y[i])
    {
      j = i;
      while ((y[j-1] < y[i]) && (j>1))
      {
        j=j-1;
        y[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
        for (m=j;m<i;m++)    y[m] = y[i];
      }
    }
  }
}





