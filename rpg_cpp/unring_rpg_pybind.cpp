/*  
    LICENCE
        
    The Software remains the property of the University Medical Center 
    Freiburg ("the University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. */

// Core functions by Elias Kellner ("unring_1D","unring_2d","Unring"),
// and Hong Hsi Lee ("unring_2d_y","unring_2d_y_2" + Partial Fourier strategy )
// Adapted to command line instruction by Ricardo Coronado-Leija 13-Feb-2023

#include <math.h>
#include <map>
#include <string.h>
#include <iostream>
#include <fstream>
#include <unistd.h>    // small 
#include <getopt.h>    // large 
#include "fftw3.h"     // fast fourier transform
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/functional.h> 
// #include <functional>

#define PI  3.1416

using namespace std;

namespace py = pybind11;

void unring_1D(fftw_complex *data,int n, int numlines,int nsh,int minW, int maxW)
{
    
    
    fftw_complex *in, *out;
    fftw_plan p,pinv;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);    
    p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    pinv = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    fftw_complex *sh = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n *(2*nsh+1));
    fftw_complex *sh2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n *(2*nsh+1));
    
    
    double nfac = 1/double(n);
    
    int *shifts = (int*) malloc(sizeof(int)*(2*nsh+1));
    shifts[0] = 0;
    for (int j = 0; j < nsh;j++)
    {
        shifts[j+1] = j+1;
        shifts[1+nsh+j] = -(j+1);
    }
    
    double TV1arr[2*nsh+1]; 
    double TV2arr[2*nsh+1]; 
    
    for (int k = 0; k < numlines; k++)
    {
     
        
        fftw_execute_dft(p,&(data[n*k]),sh);
        
        int maxn = (n%2 == 1)? (n-1)/2 : n/2 -1;
        
        for (int j = 1; j < 2*nsh+1; j++)
        {
            double phi = PI/double(n) * double(shifts[j])/double(nsh);
            fftw_complex u = {cos(phi),sin(phi)};
            fftw_complex e = {1,0};

            sh[j*n ][0] = sh[0][0];
            sh[j*n ][1] = sh[0][1];

            if (n%2 == 0)
            {
                sh[j*n + n/2][0] = 0;
                sh[j*n + n/2][1] = 0;
            }
            
            for (int l = 0; l < maxn; l++)
            {            
                
                double tmp = e[0];
                e[0] = u[0]*e[0] - u[1]*e[1];
                e[1] = tmp*u[1] + u[0]*e[1];
                
                int L ;
                L = l+1;
                sh[j*n +L][0] = (e[0]*sh[L][0] - e[1]*sh[L][1]);
                sh[j*n +L][1] = (e[0]*sh[L][1] + e[1]*sh[L][0]);
                L = n-1-l;
                sh[j*n +L][0] = (e[0]*sh[L][0] + e[1]*sh[L][1]);
                sh[j*n +L][1] = (e[0]*sh[L][1] - e[1]*sh[L][0]);                
                                
            }
        }                
        
                
        for (int j = 0; j < 2*nsh+1; j++)
        {
            fftw_execute_dft(pinv,&(sh[j*n]),&sh2[j*n]);
        }
  
        for (int j=0;j < 2*nsh+1;j++)
        {
            TV1arr[j] = 0;
            TV2arr[j] = 0;           
            const int l = 0;
            for (int t = minW; t <= maxW;t++)                
            {
                TV1arr[j] += fabs(sh2[j*n + (l-t+n)%n ][0] - sh2[j*n + (l-(t+1)+n)%n ][0]);
                TV1arr[j] += fabs(sh2[j*n + (l-t+n)%n ][1] - sh2[j*n + (l-(t+1)+n)%n ][1]);
                TV2arr[j] += fabs(sh2[j*n + (l+t+n)%n ][0] - sh2[j*n + (l+(t+1)+n)%n ][0]);
                TV2arr[j] += fabs(sh2[j*n + (l+t+n)%n ][1] - sh2[j*n + (l+(t+1)+n)%n ][1]);
            }
        }
                  

        
        
        for(int l=0; l < n; l++)
        {
            double minTV = 999999999999;
            int minidx= 0;
            for (int j=0;j < 2*nsh+1;j++)
            {                                                                    
                
//                 double TV1 = 0;
//                 double TV2 = 0;
//                 for (int t = minW; t <= maxW;t++)                
//                 {
//                     TV1 += fabs(sh2[j*n + (l-t)%n ][0] - sh2[j*n + (l-(t+1))%n ][0]);
//                     TV1 += fabs(sh2[j*n + (l-t)%n ][1] - sh2[j*n + (l-(t+1))%n ][1]);
//                     TV2 += fabs(sh2[j*n + (l+t)%n ][0] - sh2[j*n + (l+(t+1))%n ][0]);
//                     TV2 += fabs(sh2[j*n + (l+t)%n ][1] - sh2[j*n + (l+(t+1))%n ][1]);
// 
//                 }
//                   
//                 
//                 if (TV1 < minTV)
//                 {
//                     minTV = TV1;
//                     minidx = j;
//                 }
//                 if (TV2 < minTV)
//                 {
//                     minTV = TV2;
//                     minidx = j;
//                 }
                
                
                if (TV1arr[j] < minTV)
                {
                    minTV = TV1arr[j];
                    minidx = j;
                }
                if (TV2arr[j] < minTV)
                {
                    minTV = TV2arr[j];
                    minidx = j;
                }
                
                TV1arr[j] += fabs(sh2[j*n + (l-minW+1+n)%n ][0] - sh2[j*n + (l-(minW)+n)%n ][0]);
                TV1arr[j] -= fabs(sh2[j*n + (l-maxW+n)%n ][0] - sh2[j*n + (l-(maxW+1)+n)%n ][0]);
                TV2arr[j] += fabs(sh2[j*n + (l+maxW+1+n)%n ][0] - sh2[j*n + (l+(maxW+2)+n)%n ][0]);
                TV2arr[j] -= fabs(sh2[j*n + (l+minW+n)%n ][0] - sh2[j*n + (l+(minW+1)+n)%n ][0]);
                
                TV1arr[j] += fabs(sh2[j*n + (l-minW+1+n)%n ][1] - sh2[j*n + (l-(minW)+n)%n ][1]);
                TV1arr[j] -= fabs(sh2[j*n + (l-maxW+n)%n ][1] - sh2[j*n + (l-(maxW+1)+n)%n ][1]);
                TV2arr[j] += fabs(sh2[j*n + (l+maxW+1+n)%n ][1] - sh2[j*n + (l+(maxW+2)+n)%n ][1]);
                TV2arr[j] -= fabs(sh2[j*n + (l+minW+n)%n ][1] - sh2[j*n + (l+(minW+1)+n)%n ][1]);
            
            }
             
           
            double a0r = sh2[minidx*n + (l-1+n)%n ][0];
            double a1r = sh2[minidx*n + l][0];
            double a2r = sh2[minidx*n + (l+1+n)%n ][0];
            double a0i = sh2[minidx*n + (l-1+n)%n ][1];
            double a1i = sh2[minidx*n + l][1];
            double a2i = sh2[minidx*n + (l+1+n)%n ][1];
            double s = double(shifts[minidx])/nsh/2;
            
            //data[k*n + l][0] =  (a1r - 0.5*(a2r-a0r)*s + (0.5*(a2r+a0r) - a1r)*s*s)*nfac;
            //data[k*n + l][1] =  (a1i - 0.5*(a2i-a0i)*s + (0.5*(a2i+a0i) - a1i)*s*s)*nfac;

            
            if (s>0)
            {
                data[k*n + l][0] =  (a1r*(1-s) + a0r*s)*nfac;
                data[k*n + l][1] =  (a1i*(1-s) + a0i*s)*nfac;
            }
            else
            {
                s = -s;
                data[k*n + l][0] =  (a1r*(1-s) + a2r*s)*nfac;
                data[k*n + l][1] =  (a1i*(1-s) + a2i*s)*nfac;
            }
            
        }
        
       
        
    }
    
    
    
     free(shifts);
     fftw_destroy_plan(p);
     fftw_destroy_plan(pinv);
     fftw_free(in); 
     fftw_free(out);
     fftw_free(sh);
     fftw_free(sh2);
    
    
    
    
}

// Regular 2D local subvoxel-shift 
void unring_2d(fftw_complex *data1,fftw_complex *tmp2, const int *dim_sz, int nsh, int minW, int maxW)
{

    
        double eps = 0;
        fftw_complex *tmp1 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);        
        fftw_complex *data2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);    
        
        fftw_plan p,pinv,p_tr,pinv_tr;
        p = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_BACKWARD, FFTW_ESTIMATE);        
        p_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_BACKWARD, FFTW_ESTIMATE);
        double nfac = 1/double(dim_sz[0]*dim_sz[1]);
        
        for (int k = 0 ; k < dim_sz[1];k++)
           for (int j = 0 ; j < dim_sz[0];j++)
           {
                data2[j*dim_sz[1]+k][0] = data1[k*dim_sz[0]+j][0];
                data2[j*dim_sz[1]+k][1] = data1[k*dim_sz[0]+j][1];
           }
        
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);
        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
            double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0] * ck) / (ck+cj);        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1] * ck) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][0] = nfac*(tmp2[j*dim_sz[1]+k][0] * cj) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][1] = nfac*(tmp2[j*dim_sz[1]+k][1] * cj) / (ck+cj);        
            }
        }
        
        fftw_execute_dft(pinv,tmp1,data1);
        fftw_execute_dft(pinv_tr,tmp2,data2);
        
        unring_1D(data1,dim_sz[0],dim_sz[1],nsh,minW,maxW);
        unring_1D(data2,dim_sz[1],dim_sz[0],nsh,minW,maxW);
         
  
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);

        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
            // double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                // double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) ;        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) ;                                        
     //           tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) /(ck+cj);        
     //           tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) /(ck+cj);                                        
//                 tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]*ck  + tmp2[j*dim_sz[1]+k][0]*cj ) /(ck+cj);        
//                 tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]*ck  + tmp2[j*dim_sz[1]+k][1]*cj ) /(ck+cj);                                        
            }
        }
        
        fftw_execute_dft(pinv,tmp1,tmp2);
                
        fftw_free(data2);
        fftw_free(tmp1);
}

// 2D local subvoxel-shift for PF = 5/8 and = 7/8
void unring_2d_y(fftw_complex *data1,fftw_complex *tmp2, const int *dim_sz, int nsh, int minW, int maxW)
{

    
        // double eps = 0;
        fftw_complex *tmp1 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);        
        fftw_complex *data2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);    
        
        fftw_plan p,pinv,p_tr,pinv_tr;
        p = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_BACKWARD, FFTW_ESTIMATE);        
        p_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_BACKWARD, FFTW_ESTIMATE);
        double nfac = 1/double(dim_sz[0]*dim_sz[1]);
        
        for (int k = 0 ; k < dim_sz[1];k++)
           for (int j = 0 ; j < dim_sz[0];j++)
           {
                data2[j*dim_sz[1]+k][0] = data1[k*dim_sz[0]+j][0];
                data2[j*dim_sz[1]+k][1] = data1[k*dim_sz[0]+j][1];
           }
        
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);
        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
            double ck = 0;//(1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                double cj = 1;//(1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0] * ck) / (ck+cj);        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1] * ck) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][0] = nfac*(tmp2[j*dim_sz[1]+k][0] * cj) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][1] = nfac*(tmp2[j*dim_sz[1]+k][1] * cj) / (ck+cj);        
            }
        }
        
        fftw_execute_dft(pinv,tmp1,data1);
        fftw_execute_dft(pinv_tr,tmp2,data2);
        
//         unring_1D(data1,dim_sz[0],dim_sz[1],nsh,minW,maxW);
        unring_1D(data2,dim_sz[1],dim_sz[0],nsh,minW,maxW);
         
  
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);

        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
            // double ck = 0;//(1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                // double cj = 1;//(1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) ;        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) ;                                        
     //           tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) /(ck+cj);        
     //           tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) /(ck+cj);                                        
//                 tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]*ck  + tmp2[j*dim_sz[1]+k][0]*cj ) /(ck+cj);        
//                 tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]*ck  + tmp2[j*dim_sz[1]+k][1]*cj ) /(ck+cj);                                        
            }
        }
        
        fftw_execute_dft(pinv,tmp1,tmp2);
                
        fftw_free(data2);
        fftw_free(tmp1);
}

// 2D local subvoxel-shift for PF = 6/8
void unring_2d_y_2(fftw_complex *data1,fftw_complex *tmp2, const int *dim_sz, int nsh, int minW, int maxW)
{

    
        double eps = 0;
        fftw_complex *tmp1 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);       
        fftw_complex *data2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);
        
        int dim_1 = ceil(dim_sz[1]/2), dim_2 = floor(dim_sz[1]/2);
        fftw_complex *data2_1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_1);
        fftw_complex *data2_2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_2);
                
        fftw_plan p,pinv,p_tr,pinv_tr;
        p = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_BACKWARD, FFTW_ESTIMATE);        
        p_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        
        double nfac = 1/double(dim_sz[0]*dim_sz[1]);
        // double nfac_1 = 1/double(dim_sz[0]*dim_1), nfac_2 = 1/double(dim_sz[0]*dim_2);
        
        
        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
           for (int j = 0 ; j < dim_sz[0];j++)
           {
                data2[j*dim_sz[1]+k][0] = data1[k*dim_sz[0]+j][0];
                data2[j*dim_sz[1]+k][1] = data1[k*dim_sz[0]+j][1];
           }
        }
        

        
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);
        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
            double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0] * ck) / (ck+cj);        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1] * ck) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][0] = nfac*(tmp2[j*dim_sz[1]+k][0] * cj) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][1] = nfac*(tmp2[j*dim_sz[1]+k][1] * cj) / (ck+cj);        
            }
        }
        
        
        fftw_execute_dft(pinv,tmp1,data1);
        fftw_execute_dft(pinv_tr,tmp2,data2);
        
        unring_1D(data1,dim_sz[0],dim_sz[1],nsh,minW,maxW);
        unring_1D(data2,dim_sz[1],dim_sz[0],nsh,minW,maxW);
        
        
        for (int k = 0; k < dim_1; k++)
        {
            for (int j = 0; j < dim_sz[0]; j++)
            {
                data2_1[j*dim_1+k][0] = data2[j*dim_sz[1]+2*k][0];
                data2_1[j*dim_1+k][1] = data2[j*dim_sz[1]+2*k][1];;
            }
        }
        
        for (int k = 0; k < dim_2; k++)
        {
            for (int j = 0; j < dim_sz[0]; j++)
            {
                data2_2[j*dim_2+k][0] = data2[j*dim_sz[1]+2*k+1][0];
                data2_2[j*dim_2+k][1] = data2[j*dim_sz[1]+2*k+1][1];
            }
        }
        
        
        unring_1D(data2_1,dim_1,dim_sz[0],nsh,minW,maxW);
        unring_1D(data2_2,dim_2,dim_sz[0],nsh,minW,maxW);
        
        for (int k = 0; k < dim_1; k++)
        {
            for (int j = 0; j< dim_sz[0]; j++)
            {
                data2[j*dim_sz[1]+2*k][0] = data2_1[j*dim_1+k][0];
                data2[j*dim_sz[1]+2*k][1] = data2_1[j*dim_1+k][1];
            }
        }
        
        for (int k = 0; k < dim_2; k++)
        {
            for (int j = 0; j< dim_sz[0]; j++)
            {
                data2[j*dim_sz[1]+2*k+1][0] = data2_2[j*dim_2+k][0];
                data2[j*dim_sz[1]+2*k+1][1] = data2_2[j*dim_2+k][1];
            }
        }
        
        
  
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);

        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
            // double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                // double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) ;        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) ;                                        
     //           tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) /(ck+cj);        
     //           tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) /(ck+cj);                                        
//                 tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]*ck  + tmp2[j*dim_sz[1]+k][0]*cj ) /(ck+cj);        
//                 tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]*ck  + tmp2[j*dim_sz[1]+k][1]*cj ) /(ck+cj);                                        
            }
        }
        
        fftw_execute_dft(pinv,tmp1,tmp2);
                
        fftw_free(data2);
        fftw_free(tmp1);
        fftw_free(data2_1);
        fftw_free(data2_2);

}


// Applying unrining to 2D images ( dimensions > 2 are not used, so that part was removed )
void Unring(double *data, double *data_i, double *res, double *res_i,  int *dim_sz, 
unsigned int numdim, unsigned int pfo, unsigned int minW, unsigned int maxW, unsigned int nsh){

                
        fftw_complex *data_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);    
        fftw_complex *res_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);    
        if (data_i == 0)
        {
            // plhs[0] = mxCreateNumericArray(numdim,dim_sz,mxGetClassID(Img),mxREAL);
            for (int k = 0 ; k < dim_sz[1];k++)
               for (int j = 0 ; j < dim_sz[0];j++)
               {
                    data_complex[k*dim_sz[0]+j][0] = data[k*dim_sz[0]+j];
                    data_complex[k*dim_sz[0]+j][1] = 0;
               }                        
        }
        else
        {
            // plhs[0] = mxCreateNumericArray(numdim,dim_sz,mxGetClassID(Img),mxCOMPLEX);
            for (int k = 0 ; k < dim_sz[1];k++)
               for (int j = 0 ; j < dim_sz[0];j++)
               {
                    data_complex[k*dim_sz[0]+j][0] = data[k*dim_sz[0]+j];
                    data_complex[k*dim_sz[0]+j][1] = data_i[k*dim_sz[0]+j];
               }                        
        }
        
        if(pfo == 1){
        unring_2d_y(data_complex,res_complex,dim_sz,nsh,minW,maxW);
        }
        else if(pfo == 2){
        unring_2d_y_2(data_complex,res_complex,dim_sz,nsh,minW,maxW);
        }
        else if(pfo == 3){
        unring_2d_y(data_complex,res_complex,dim_sz,nsh,minW,maxW);
        }
        else{
        unring_2d(data_complex,res_complex,dim_sz,nsh,minW,maxW);    
        }
        
        
        // double *res =  (double*) mxGetData(plhs[0]);
        // double *res_i =  (double*) mxGetImagData(plhs[0]);
        
        if (res_i == 0)
        {
            for (int k = 0 ; k < dim_sz[1];k++)
               for (int j = 0 ; j < dim_sz[0];j++)
                   res[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][0];
        } 
        else
        {
            for (int k = 0 ; k < dim_sz[1];k++)
               for (int j = 0 ; j < dim_sz[0];j++)
               {
                   res[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][0];
                   res_i[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][1];
               }
        } 
                        
        fftw_free(data_complex);
        fftw_free(res_complex);

}

// ===== start useful functions ===== //
template <typename T>
T **Matrix(unsigned int nr, unsigned int nc){
   T **Matrix;
   unsigned int i;
   Matrix = new T *[nr];
   if(Matrix != NULL){
      for(i = 0; i < nr; i++){
         Matrix[i] = new T[nc];
      }
   } else {
      std::cout << "ERROR: matrix could not be created" << std::endl; 
      exit(1);
   }
   return Matrix;
}

template
double **Matrix(unsigned int nr, unsigned int nc);


template <typename T>
void FreeMatrix(T **M, int nr, int nc){
   for(int i = 0; i < nr; i++) delete[] M[i];
   delete[] M;
}


template <typename T>
void Reshape(T* volume, T** slices, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int nvol, bool pfdimf, bool toSlices) {
    unsigned int nxy = nx * ny;
    unsigned int nxyz = nxy * nz;
    unsigned int idx_volume, idx_slices;

    if (toSlices) {
        for (unsigned int l = 0; l < nvol; ++l) {
            for (unsigned int k = 0; k < nz; ++k) {
                for (unsigned int j = 0; j < ny; ++j) {
                    for (unsigned int i = 0; i < nx; ++i) {
                        if (pfdimf) {
                            // Reordering from [nvol, nslice, ny, nx] to [nx, ny, nslice, nvol]
                            idx_volume = l * nxyz + k * nxy + j * nx + i;
                            idx_slices = j * nx + i;
                        } else {
                            // Reordering from [nvol, nslice, ny, nx] to [nx, ny, nslice, nvol]
                            idx_volume = l * nxyz + k * nxy + j * nx + i;
                            idx_slices = i * ny + j;
                        }
                        slices[l * nz + k][idx_slices] = volume[idx_volume];
                    }
                }
            }
        }
    } else {
        for (unsigned int l = 0; l < nvol; ++l) {
            for (unsigned int k = 0; k < nz; ++k) {
                for (unsigned int j = 0; j < ny; ++j) {
                    for (unsigned int i = 0; i < nx; ++i) {
                        if (pfdimf) {
                            // Reordering from [nx, ny, nslice, nvol] back to [nvol, nslice, ny, nx]
                            idx_volume = l * nxyz + k * nxy + j * nx + i;
                            idx_slices = j * nx + i;
                        } else {
                            // Reordering from [nx, ny, nslice, nvol] back to [nvol, nslice, ny, nx]
                            idx_volume = l * nxyz + k * nxy + j * nx + i;
                            idx_slices = i * ny + j;
                        }
                        volume[idx_volume] = slices[l * nz + k][idx_slices];
                    }
                }
            }
        }
    }
}
template <typename T>
void Reshape(T* volume, T** slices, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int nvol, bool pfdimf, bool toSlices);


// MATLABs NN
template <typename T>
void NN(T *I, T *O, unsigned int ini, unsigned int inx, unsigned int iny, unsigned int inio, unsigned int onx, unsigned int ony){

unsigned int x, y;
unsigned int src_x, src_y;
double ux, uy, scale_x, scale_y;

scale_x = onx/(double)inx;
scale_y = ony/(double)iny;
// cout << scale_x << " " << scale_y << endl;

// imresize.m -> contributions.m
// Input-space coordinates. Calculate the inverse mapping such that 0.5
// in output space maps to 0.5 in input space, and 0.5+scale in output
// space maps to 1.5 in input space.
for(x = 0; x < onx; x++){
ux    = (x+1)/scale_x + 0.5*(1.0 - 1.0/scale_x); // here I added x -> x+1        
src_x = floor( ux - 0.5 );                       // here I removed +1
for(y = 0; y < ony; y++){
uy    = (y+1)/scale_y + 0.5*(1.0 - 1.0/scale_y);    
src_y = floor( uy - 0.5 );    
O[inio + y*onx + x] = I[ini + src_y*inx + src_x];
// cout << ux << " " << uy << " " << src_x << " " << src_y << endl;
} // x    
} // y
// cout << endl;
// cout << endl;

} // NN
template
void NN(double *I, double *O, unsigned int ini , unsigned int inr, unsigned int inc, unsigned int inio, unsigned int onr, unsigned int onc);


template <typename T>
void UnringScale(T *I, T *Ii, unsigned int nx, unsigned int ny, unsigned int scale, double yfact, unsigned int minW, unsigned int maxW, unsigned int nsh){

T *O;
T *Oi;
bool pf7_8 = fabs(yfact - 1) > 1e-6; // here yfact should be 1 (5/8) or 3 (7/8)
bool fimag = (Ii != nullptr);
//printf("UnringScale: yfact = %f, pf7_8 =  %d, fimag =  %d\n",yfact,pf7_8,fimag);
unsigned int i, nyy;
if(pf7_8){ // 7/8
nyy = (int)round(yfact*ny);    
O   = new T[nx*nyy];    
NN(I,O,0,nx,ny,0,nx,nyy);
if(fimag){
Oi  = new T[nx*nyy];    
NN(Ii,Oi,0,nx,ny,0,nx,nyy);
} // fimag
} // 7/8 => NN
else{ // 5/8
nyy = ny;    
O   = new T[nx*ny];    
for(i = 0; i < nx*ny; i++){
O[i] = I[i];
} // i  
if(fimag){
Oi  = new T[nx*ny];    
for(i = 0; i < nx*ny; i++){
Oi[i] = Ii[i];
} // i      
} // fimag
} // 5/8
// WriteImage(O,nx,nyy,(char*)"image.pgm");

unsigned int j, x, y, nys, pfo, ndim;
unsigned int nys_[scale];
int *dim_sz = new int[2];
nys  = ceil(nyy/(double)scale);
pfo  = 0; // regular sushi
ndim = 2;

double **Is = Matrix<double>(scale,nx*nys);
double **Os = Matrix<double>(scale,nx*nys);
double **Isi, **Osi;
if(fimag){
Isi = Matrix<double>(scale,nx*nys);
Osi = Matrix<double>(scale,nx*nys);
}// fimag
else{
Isi = 0;
Osi = 0;    
}

for(i = 0; i < scale; i++){    
j = 0;    
for(y = i; y < nyy; y+=scale){  
for(x = 0; x < nx;  x++){  
Is[i][j*nx+x] = O[y*nx+x];
if(fimag){
Isi[i][j*nx+x] = Oi[y*nx+x];
} // fimag
} // x
j++;
} // y

nys_[i]   = j; 
dim_sz[0] = nx; 
dim_sz[1] = nys_[i];
if(fimag){
Unring(Is[i],Isi[i],Os[i],Osi[i],dim_sz,ndim,pfo,minW,maxW,nsh);    
}
else{
Unring(Is[i],0,Os[i],0,dim_sz,ndim,pfo,minW,maxW,nsh);    
}

j = 0;    
for(y = i; y < nyy; y+=scale){  
for(x = 0; x < nx;  x++){  
O[y*nx+x] = Os[i][j*nx+x];
if(fimag){
Oi[y*nx+x] = Osi[i][j*nx+x];
} // fimag
} // x
j++;
} // y
} // i

// cout << nys << " " << nys_[0] << " " << nys_[1] << " " << nys_[2] << " " << nys_[3] << endl;

if(pf7_8){   
NN(O,I,0,nx,nyy,0,nx,ny);
if(fimag){
NN(Oi,Ii,0,nx,nyy,0,nx,ny);    
} // fimag
} // NN
else{
for(i = 0; i < nx*ny; i++){
if(fimag){    
Ii[i] = Oi[i];
} // fimag
} // i   
}

delete[] O;
delete[] dim_sz;
FreeMatrix(Is,scale,nx*nys);
FreeMatrix(Os,scale,nx*nys);
if(fimag){
delete[] Oi;    
FreeMatrix(Isi,scale,nx*nys);
FreeMatrix(Osi,scale,nx*nys);
} // fimag
} // UnringScale
template
void UnringScale(double *I, double *Ii, unsigned int nx, unsigned int ny, unsigned int scale, double yfact, unsigned int minW, unsigned int maxW, unsigned int nsh);

// ===== end useful functions ===== //

py::tuple unring(py::array_t<double> data, py::array_t<double> phase = py::array_t<double>(), 
                  unsigned int minW = 1, unsigned int maxW = 3, unsigned int nsh = 20, 
                  double pfv = 6.0 / 8.0, bool pfdimf = true, bool phase_flag = false) {
    py::buffer_info data_info = data.request();
    py::buffer_info phase_info = phase_flag ? phase.request() : py::buffer_info();

    double* volumein = static_cast<double*>(data_info.ptr);
    double* phasein = phase_flag ? static_cast<double*>(phase_info.ptr) : nullptr;

    unsigned int nx = data_info.shape[3];
    unsigned int ny = data_info.shape[2];
    unsigned int nz = data_info.shape[1];
    unsigned int ndwi = data_info.shape[0];
    unsigned int nelem = nx * ny * nz * ndwi;

    std::vector<double> volumeout(nelem, 0.0);
    std::vector<double> phaseout(phase_flag ? nelem : 0, 0.0);

    // Absolute and rescaling
    for (unsigned int i = 0; i < nelem; ++i) {
        volumein[i] = std::abs(volumein[i]);
    }

    if (phase_flag) {
        for (unsigned int i = 0; i < nelem; ++i) {
            // phasein[i] = std::abs(phasein[i]);
            // Check phase in radians
            if (std::abs(phasein[i]) > 3.1416) {
                throw std::runtime_error("Error: Phase contains values outside the range [-pi, pi].");
            }
        }

        // Compute complex (replace magnitude/phase)
        for (unsigned int i = 0; i < nelem; ++i) {
            double mag = volumein[i];
            double pha = phasein[i];
            volumein[i] = mag * std::cos(pha);
            phasein[i] = mag * std::sin(pha);
        }
    }

    // Create matrices for reshaping
    double** slicesin = Matrix<double>(nz * ndwi, nx * ny);
    double** slicesout = Matrix<double>(nz * ndwi, nx * ny);
    Reshape<double>(volumein, slicesin, nx, ny, nz, ndwi, pfdimf, true);

    double** slicesin_i = nullptr;
    double** slicesout_i = nullptr;
    if (phase_flag) {
        slicesin_i = Matrix<double>(nz * ndwi, nx * ny);
        slicesout_i = Matrix<double>(nz * ndwi, nx * ny);
        Reshape<double>(phasein, slicesin_i, nx, ny, nz, ndwi, pfdimf, true);
    }

    // Set dimensions
    int dim_sz[4] = { static_cast<int>(nx), static_cast<int>(ny), static_cast<int>(nz), static_cast<int>(ndwi) };
    if (pfdimf) { // y
        std::swap(dim_sz[0], dim_sz[1]);
    }

    // Determine pfo based on pfv
    int pfo;
    if (std::abs(pfv - 5.0 / 8.0) < 1e-6) {
        pfo = 1;
    } else if (std::abs(pfv - 6.0 / 8.0) < 1e-6) {
        pfo = 2;
    } else if (std::abs(pfv - 7.0 / 8.0) < 1e-6) {
        pfo = 3;
    } else if (std::abs(pfv - 1.0) < 1e-6) {
        pfo = 0;
    } else {
        throw std::runtime_error("Error: this command only supports Full Fourier and PF = 5/8, 6/8, and 7/8.");
    }

    std::cout << "minW: " << minW << std::endl;
    std::cout << "maxW: " << maxW << std::endl;
    std::cout << "nsh: " << nsh << std::endl;
    std::cout << "partial fourier fraction: " << pfv << std::endl;
    std::cout << "partial fourier dimension: " << pfdimf << std::endl;
    std::cout << "phase flag: " << phase_flag << std::endl;
    std::cout << "number of slices: " << nz * ndwi << std::endl;

    // Apply unringing
    int scale = 4;
    for (unsigned int i = 0; i < nz * ndwi; ++i) {
        if (PyErr_CheckSignals() != 0) {
            throw py::error_already_set();  // Raise the Python exception if a signal was received
        }

        if (pfo == 1 || pfo == 3) { // Cases where ifact can be equal to pfo
            if (pfdimf) { // y
                if (phase_flag) { // complex
                    UnringScale(slicesin[i], slicesin_i[i], nx, ny, scale, pfo, minW, maxW, nsh);
                } else { // real
                    UnringScale(slicesin[i], static_cast<double*>(nullptr), nx, ny, scale, pfo, minW, maxW, nsh);
                }
            } else { // x
                if (phase_flag) { // complex
                    UnringScale(slicesin[i], slicesin_i[i], ny, nx, scale, pfo, minW, maxW, nsh);
                } else { // real
                    UnringScale(slicesin[i], static_cast<double*>(nullptr), ny, nx, scale, pfo, minW, maxW, nsh);
                }
            }
        }
        if (phase_flag) { // complex
            Unring(slicesin[i], slicesin_i[i], slicesout[i], slicesout_i[i], dim_sz, 2, pfo, minW, maxW, nsh);
        } else { // real
            Unring(slicesin[i], 0, slicesout[i], 0, dim_sz, 2, pfo, minW, maxW, nsh);
        }
    }

    // Put slices into volume (real)
    Reshape<double>(volumeout.data(), slicesout, nx, ny, nz, ndwi, pfdimf, false);
    FreeMatrix(slicesin, nz * ndwi, nx * ny);
    FreeMatrix(slicesout, nz * ndwi, nx * ny);

    if (phase_flag) { // (imag)
        Reshape<double>(phaseout.data(), slicesout_i, nx, ny, nz, ndwi, pfdimf, false);
        FreeMatrix(slicesin_i, nz * ndwi, nx * ny);
        FreeMatrix(slicesout_i, nz * ndwi, nx * ny);

        // Compute magnitude/phase (replace real/imag)
        for (unsigned int i = 0; i < nelem; ++i) {
            double re = volumeout[i];
            double im = phaseout[i];
            volumeout[i] = sqrt(re * re + im * im);
            phaseout[i] = atan2(im, re);
        }
    }

    return py::make_tuple(
        py::array_t<double>(data_info.shape, volumeout.data()),
        phase_flag ? py::array_t<double>(phase_info.shape, phaseout.data()) : py::array_t<double>()
    );
}


PYBIND11_MODULE(rpg, m) {
    m.def("unring", &unring, "Process data with given parameters",
          py::arg("data"), py::arg("phase") = py::array_t<double>(),
          py::arg("minW") = 1, py::arg("maxW") = 3, py::arg("nsh") = 20,
          py::arg("pfv") = 6.0 / 8.0, py::arg("pfdimf") = true, py::arg("phase_flag") = false);
}