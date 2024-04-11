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
#include <nifti1_io.h> // medical images
#include "fftw3.h"     // fast fourier transform

#define PI  3.1416

using namespace std;

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
            double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
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
            double ck = 0;//(1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                double cj = 1;//(1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
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
        double nfac_1 = 1/double(dim_sz[0]*dim_1), nfac_2 = 1/double(dim_sz[0]*dim_2);
        
        
        
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
            double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
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

// ===== start nifti functionality ===== //
template <typename Tin, typename Tout>
void PtrVoid2Double(void *nii, Tout *data, int N){
   Tin *d = reinterpret_cast<Tin*> (nii); 
   for(int i = 0; i < N; i++) data[i] = (Tout)d[i];
   }

template <typename T>
void GetNiftiData(void *nii, T *data, int N, short int datatype){
	
switch(datatype){
   case NIFTI_TYPE_INT8:    PtrVoid2Double<char,T>(nii,data,N);           break; // 256
   case NIFTI_TYPE_UINT8:   PtrVoid2Double<unsigned char,T>(nii,data,N);  break; // 2
   case NIFTI_TYPE_INT16:   PtrVoid2Double<short,T>(nii,data,N);          break; // 4
   case NIFTI_TYPE_UINT16:  PtrVoid2Double<unsigned short,T>(nii,data,N); break; // 512
   case NIFTI_TYPE_INT32:   PtrVoid2Double<int,T>(nii,data,N);            break; // 8
   case NIFTI_TYPE_UINT32:  PtrVoid2Double<unsigned int,T>(nii,data,N);   break; // 768
   case NIFTI_TYPE_FLOAT32: PtrVoid2Double<float,T>(nii,data,N);          break; // 16
   case NIFTI_TYPE_FLOAT64: PtrVoid2Double<double,T>(nii,data,N);         break; // 64
   default: cout << "GetNiftiDataDouble() ERROR: datatype not valid " << datatype <<endl;        
   }
   
}

template
void GetNiftiData(void *nii, unsigned short *data, int N, short int datatype);
template
void GetNiftiData(void *nii, float *data, int N, short int datatype);
template
void GetNiftiData(void *nii, double *data, int N, short int datatype);
// ===== end nifti functionality ===== //

// ===== start string operations ===== //

// convert any uppercase chars to lowercase
int MakeLowerCase(char *str){
   int c;

   if(!str || !*str) return 0;

   for(c = 0; c < strlen(str); c++)
      if( isupper(str[c]) ) str[c] = tolower(str[c]);

   return 0;
}

// Run strcmp against of list of strings
// return index of equality, if found
// else return -1 
int CompareStringList(const char *str, char **strlist, int len){
   int c;
   if(len <= 0 || !str || !strlist) return -1;
   for(c = 0; c < len; c++)
      if(strlist[c] && !strcmp(str, strlist[c])) return c;
   return -1;
   }

// Duplicate the given string (alloc length+1)
// return allocated pointer (or NULL on failure)
char *StringDuplicate(const char *str){
   char *dup;

   // allow calls passing NULL
   if(!str) return NULL;        
   
   dup = (char *)malloc(strlen(str) + 1);
   // check for failure
   if(dup) strcpy(dup,str);
   else    fprintf(stderr,"StringDuplicate() ERROR: failed to alloc %u bytes\n", 
                   (unsigned int)strlen(str)+1);

   return dup;
   }

// Check the end of the filename for a valid file extension
//    Valid extensions are currently .nii, .hdr, .img, .nia,
//    or any of them followed by .gz.  Note that '.' is part of
//    the extension.
//
//    Uppercase extensions are also valid, but not mixed case.
//
//    \return a pointer to the extension (within the filename), or NULL
char *FindFileExtension(const char *name){
   char *ext, extcopy[8];
   int    len;
   char   extnii[8] = ".nii";   // modifiable, for possible uppercase
   char   exthdr[8] = ".hdr";   
   char   extimg[8] = ".img";
   char   extnia[8] = ".nia";
   char   exttxt[8] = ".txt";
   char   extdat[8] = ".dat";
   char *elist[6]  = { NULL, NULL, NULL, NULL, NULL, NULL};

   // stupid compiler... 
   elist[0] = extnii; elist[1] = exthdr; elist[2] = extimg; elist[3] = extnia;
   elist[4] = exttxt; elist[5] = extdat;

   if (!name) return NULL;

   len = (int)strlen(name);
   if (len < 4) return NULL;

   ext = (char *)name + len - 4;

   // make manipulation copy, and possibly convert to lowercase
   strcpy(extcopy,ext); 
   MakeLowerCase(extcopy);

   // if it look like a basic extensionreturn it
   if(CompareStringList(extcopy, elist, 6) >= 0){
      return ext;
      }

   return NULL;
   }

//Duplicate the filename, while clearing any extension
//This allocates memory for basename which should eventually be freed.
char *MakeBaseName(const char* fname){
   char *basename, *ext;

   basename = StringDuplicate(fname);

   ext = FindFileExtension(basename);
   if ( ext ) *ext = '\0';  // clear out extension
   
   return basename;  // in either case
   }

// ===== end string operations ===== //

// ===== start useful functions ===== //
template <typename T>
T **Matrix(unsigned int nr, unsigned int nc){
   T **Matrix;
   int i;
   Matrix = new T *[nr];
   if(Matrix != NULL){
      for(i = 0; i < nr; i++){
         Matrix[i] = new T[nc];
         }// i
      }// if
   else{
      cout << "ERROR: matrix could not be created" << endl; exit(1);
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
template
void FreeMatrix(double **M, int nr, int nc);

template <typename T>
void Reshape(T *I, T **O, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int ndwi, bool t, bool dir){
unsigned x, y, z, d, idxim, idxnii;
for(d = 0; d < ndwi; d++){
for(z = 0; z < nz;   z++){
for(y = 0; y < ny;   y++){
for(x = 0; x < nx;   x++){
if(t){ // dim y
idxim  = y*nx + x;    
}   
else{ // dim x
idxim  = x*ny + y;    
} 
idxnii = d*(nx*ny*nz) + z*(nx*ny) + y*nx + x;
if(dir){ // before
O[d*nz+z][idxim] = I[idxnii];
}
else{ // after
I[idxnii] = O[d*nz+z][idxim];
} 
} // x    
} // y   
} // z    
} // d

} // Reshape 
template
void Reshape(double *I, double **O, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int ndwi, bool t, bool dir);

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
bool fimag = (Ii != 0);
printf("UnringScale: yfact = %f, pf7_8 =  %d, fimag =  %d\n",yfact,pf7_8,fimag);
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

int show_help(){
printf(
// === Non Optional Arguments === //
"rpg [ options ] input output\n"
"\n"
"\tinput\n"
"\t\tname of the input volume file .nii(.gz)\n"
"\toutput\n"
"\t\tname of the output volume file .nii(.gz)\n"    
"\n"
"Tool for removal of the Gibbs ringing artefact with or without Partial Fourier in Magnetic Resonance Images.\n"
"\n"
// === Optional Arguments === //
"Options:\n"
"\n"
"\t -phase file\n"
"\t\tname for the input phase volume file .nii(.gz), if available. Should be in the range [-pi,pi] radians.\n" 
"\t -pf option\n"
"\t\tpartial Fourier factor, support only for PF = 5/8, 6/8, 7/8, and 1 (no PF). Default: 6/8\n" 
"\t -dim option\n"
"\t\tpartial Fourier dimension (x, 0, i) or (y, 1, j). Do not matter for PF = 1. Default: y\n"      
"\t -minw value\n"
"\t\tleft border of window used for TV computation: Default: 1.\n" 
"\t -maxw value\n"
"\t\tright border of window used for TV computation: Default: 3.\n" 
"\t -nsh value\n"
"\t\tdiscretization of subpixel spaceing: Default: 20.\n"  
"\t -help\n"
"\t\tshow this help\n"
"\n"
"References:\n"
"\tPartial Fourier (RPG):\n"
"\t\tLee, et al., Magnetic Resonance in Medicine, 2021 (doi.org/10.1002/mrm.28830)\n"
"\tFull Fourier (Local Subvoxel-Shifts):\n"
"\t\tKellner, et al., Magnetic Resonance in Medicine, 2015 (doi.org/10.1002/mrm.26054)\n"
"\n");
return 0;
}    

int main(int argc, char** argv){

// .......................................................................... //    
// ########################################################################## //
// .......................................................................... //
   
int option_index; // getopt_long_only stores the option index here.
int goloval;      // getopt_long_only returns the option value here

bool help_flag      = false;
bool phase_flag     = false;  
// Parameters for local subvoxel shifts
unsigned int minW   = 1;
unsigned int maxW   = 3;
unsigned int nsh    = 20;
// Parameters PRG
double       pfv    = 6.0/8.0; // partial fourier value
unsigned int pfo    = 2;       // partial fourier option (1=5/8, 2=6/8, 3=7/8, anyother=1) 
bool         pfdimf = true;    // partal fourier dimension false = x(or 0), true = y(or 1)
char namePhase[500];

// options
struct option long_options[] = {
  // These options don’t set a flag. We distinguish them by their indices.
  {"minw", required_argument, 0, 'w'},
  {"maxw", required_argument, 0, 'W'},      
  {"nsh",  required_argument, 0, 's'},
  {"pf",   required_argument, 0, 'p'},
  {"dim",  required_argument, 0, 'd'},
  {"phase",required_argument, 0, 'a'},
  {"help",       no_argument, 0, 'h'},
  {0, 0, 0, 0}
  };

while(true){
  // Detecting the next option
  goloval = getopt_long_only (argc, argv, "w:W:s:p:d:a:h", long_options, &option_index);
  
  // Detect the end of the options and break the while. 
  if (goloval == -1) break;

  switch (goloval){
      case 'w':
         minW  = atoi(optarg);
         break;   
      case 'W':
         maxW  = atoi(optarg);
         break; 
      case 's':
         nsh   = atoi(optarg);
         break;          
      case 'p':
         if(strcmp(optarg,"5/8") == 0){
         pfv = 5.0/8.0; 
         }
         else if(strcmp(optarg,"6/8") == 0){
         pfv = 6.0/8.0; 
         }
         else if(strcmp(optarg,"7/8") == 0){
         pfv = 7.0/8.0; 
         }
         else if(strcmp(optarg,"8/8") == 0){
         pfv = 8.0/8.0; 
         }
         else{
         pfv = atof(optarg);
         }
         // cout << optarg << " " << pfv << " " << strcmp(optarg,"5/8") <<  endl;
         break; 
      case 'd':
         if(strcmp(optarg,"0") == 0 || strcmp(optarg,"x") == 0 || strcmp(optarg,"X") == 0 || strcmp(optarg,"i") == 0 || strcmp(optarg,"I") == 0){
         pfdimf = false;  
         } 
         else if(strcmp(optarg,"1") == 0 || strcmp(optarg,"y") == 0 || strcmp(optarg,"Y") == 0 || strcmp(optarg,"J") == 0 || strcmp(optarg,"J") == 0){
         pfdimf = true;  
         }
         else{
         cout << "Error: not valid value for option -dim" << endl;   
         exit(1);
         }
         // cout << optarg << " " << pfdimf << endl;
         break; 
      case 'a':
         phase_flag = true;
         strcpy(namePhase,optarg); 
         break;         
      case 'h':
         help_flag = true;
         break;
      case '?':
         //getopt_long already printed an error message.
         break;
      default:
         abort();
      } // switch
  } // while (1)

// Showing the help
if(help_flag || argc == 1){
  show_help();
  exit(1);
  }

// ======================================================================= //	  
// Remaining command line arguments (not options).
if(argc - optind != 2){
  printf("Expecting 2 arguments %d provided. Type rpg -help for help.\n",(argc - optind));
  return 1;
  }

// --- Non-Optional arguments --- //
// Reading input parameters  
char nameDataIn[300], nameDataOut[300];

sprintf(nameDataIn ,"%s" ,argv[optind++]);
sprintf(nameDataOut,"%s" ,argv[optind++]);

// ======================================================================= //

// PF option selected (int)
if( fabs(pfv - 5.0/8.0) < 1e-6){
pfo = 1;
}
else if(fabs(pfv - 6.0/8.0) < 1e-6){
pfo = 2; 
}
else if(fabs(pfv - 7.0/8.0) < 1e-6){
pfo = 3; 
}
else if(fabs(pfv - 1.0) < 1e-6){
pfo = 0; 
}
else{
cout << "Error: this command only supports Full Fourier and PF = 5/8, 6/8, and 7/8." << endl; 
exit(1);
}

// if phase is provided create a basename from output filename, as there will be two outputs
char *nameDataOutNew, *extOut;
char nameDataOutMag[300], nameDataOutPha[300];
if(phase_flag){
extOut = FindFileExtension(nameDataOut);
nameDataOutNew = MakeBaseName(nameDataOut);    
sprintf(nameDataOutMag,"%s_mag%s",nameDataOutNew,extOut);
sprintf(nameDataOutPha,"%s_phase%s",nameDataOutNew,extOut);
}

// print options
printf("RPG\n");
if(phase_flag){
printf("Inputs:\n");    
printf("Magnitude: %s\n",nameDataIn);  
printf("Phase: %s\n",namePhase);    
}
else{
printf("Input: %s\n",nameDataIn);
}
if(phase_flag){
printf("Outputs:\nMag=%s\nPhase=%s\n",nameDataOutMag,nameDataOutPha);
}
else{
printf("Output: %s\n",nameDataOut);
}
printf("Options:\n");
printf("pfdim = %s, pf = %.3f(option %d)\n",(pfdimf ? "y" : "x"),pfv,pfo);
printf("minw = %d, maxw =  %d, nsh =  %d\n",minW,maxW,nsh);

// .......................................................................... //    
// ########################################################################## //
// .......................................................................... //

bool ifact;
unsigned int ndim, nx, ny, nz, ndwi, nvox2, nvox3, nelem;
unsigned int i, j, x, y, z, c, d, idx, idxnii, scale;

// read input volume
nifti_image    *imin = nifti_image_read(nameDataIn,true);
nifti_1_header spec  = nifti_convert_nim2nhdr(imin);
// spec.datatype        = NIFTI_TYPE_FLOAT64;
// spec.bitpix          = 64;
spec.datatype        = NIFTI_TYPE_FLOAT32;
spec.bitpix          = 32;
spec.scl_slope       = 1.0; 
spec.scl_inter       = 0.0; 

ndim  = imin->dim[0];
nx    = imin->dim[1];
ny    = imin->dim[2];
nz    = imin->dim[3];
ndwi  = imin->dim[4];
nvox2 = nx*ny;
nvox3 = nx*ny*nz;
nelem = nx*ny*nz*ndwi;

printf("ndim = %d : nx = %d, ny = %d, nz = %d, ndwi = %d\n",ndim,nx,ny,nz,ndwi);

// get data
double *volumein  = new double[nelem];
double *volumeout = new double[nelem];
GetNiftiData<double>(imin->data,volumein,nelem,imin->datatype); 
// absolute and rescaling
for(i = 0; i < nelem; i++){
if(fabs(imin->scl_slope) > 1e-12){    
volumein[i] = fabs((imin->scl_slope)*volumein[i]+(imin->scl_inter));    
}
else{
volumein[i] = fabs(volumein[i]+(imin->scl_inter));    
}
} // i
nifti_image_free(imin);

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// if phase is provided set complex data (real/imag)

double mag_, pha_, re_, im_;
bool flag_norad = false;
double *phasein;
double *phaseout;

if(phase_flag){
// read phase volume   
nifti_image *pmin = nifti_image_read(namePhase,true);
// check same size
if(ndim!=pmin->dim[0] || nx!=pmin->dim[1] || ny!=pmin->dim[2] || nz!=pmin->dim[3] || ndwi!=pmin->dim[4]){
printf("RPG ERROR: Phase should have the same dimensions as input");
exit(1);    
}

// get data
phasein  = new double[nelem];
phaseout = new double[nelem];
GetNiftiData<double>(pmin->data,phasein,nelem,pmin->datatype); 
// absolute and rescaling
for(i = 0; i < nelem; i++){
if(fabs(pmin->scl_slope) > 1e-12){    
phasein[i] = fabs((pmin->scl_slope)*phasein[i]+(pmin->scl_inter));    
}
else{
phasein[i] = fabs(phasein[i]+(pmin->scl_inter));    
}
// check phase in radians
if(fabs(phasein[i]) > 3.1416){
flag_norad = true;
} // if
} // i
if(flag_norad){
printf("RPG WARNING: Phase contains values outside the range [-pi,pi];");
} //
nifti_image_free(pmin);

// compute complex (replace magnitude/phase)
for(i = 0; i < nelem; i++){
mag_ = volumein[i];
pha_ = phasein[i];
re_  = mag_*cos(pha_);
im_  = mag_*sin(pha_);
volumein[i] = re_;
phasein[i]  = im_;
} // i

} // phase 

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

// set dimensions
int *dim_sz = new int[4];
dim_sz[2] = nz;
dim_sz[3] = ndwi;
if(pfdimf){ // dim y
dim_sz[0] = nx;
dim_sz[1] = ny;
}
else{ // dim x
dim_sz[0] = ny;
dim_sz[1] = nx;    
}

// divide volume (real) into slices
double **slicesin  = Matrix<double>(nz*ndwi,nx*ny);
double **slicesout = Matrix<double>(nz*ndwi,nx*ny);
Reshape<double>(volumein,slicesin,nx,ny,nz,ndwi,pfdimf,1);
delete[] volumein;

double **slicesin_i;
double **slicesout_i;
if(phase_flag){
// divide volume (imag) into slices    
slicesin_i  = Matrix<double>(nz*ndwi,nx*ny);
slicesout_i = Matrix<double>(nz*ndwi,nx*ny);
Reshape<double>(phasein,slicesin_i,nx,ny,nz,ndwi,pfdimf,1);
delete[] phasein;
} // phase
else{
slicesin_i = 0;
slicesout_i = 0;    
}

// apply unringing
scale = 4;
ndim  = 2;
for(i = 0; i < nz*ndwi; i++){
if(pfo == 1 || pfo == 3){ // just for these cases ifact can be equal to pfo
if(pfdimf){ // y    
if(phase_flag) // complex   
UnringScale<double>(slicesin[i],slicesin_i[i],nx,ny,scale,pfo,minW,maxW,nsh);
else // real
UnringScale<double>(slicesin[i],0,nx,ny,scale,pfo,minW,maxW,nsh);    
}
else{ // x
if(phase_flag) // complex       
UnringScale<double>(slicesin[i],slicesin_i[i],ny,nx,scale,pfo,minW,maxW,nsh);    
else // real
UnringScale<double>(slicesin[i],0,ny,nx,scale,pfo,minW,maxW,nsh);    
} // dim
} // pfo
if(phase_flag) // complex
Unring(slicesin[i],slicesin_i[i],slicesout[i],slicesout_i[i],dim_sz,ndim,pfo,minW,maxW,nsh);
else // real
Unring(slicesin[i],0,slicesout[i],0,dim_sz,ndim,pfo,minW,maxW,nsh);    
// slicesout = slicesin;
} // i

// put slices into volume (real)
Reshape<double>(volumeout,slicesout  ,nx,ny,nz,ndwi,pfdimf,0);
FreeMatrix(slicesin ,nz*ndwi,nx*ny);
FreeMatrix(slicesout,nz*ndwi,nx*ny);
if(phase_flag){ // (imag)
Reshape<double>(phaseout ,slicesout_i,nx,ny,nz,ndwi,pfdimf,0);
FreeMatrix(slicesin_i ,nz*ndwi,nx*ny);
FreeMatrix(slicesout_i,nz*ndwi,nx*ny);

// compute magnitude/phase (replace real/imag)
for(i = 0; i < nelem; i++){
re_  = volumeout[i];
im_  = phaseout[i];
mag_ = sqrt(re_*re_ + im_*im_);
pha_ = atan2(im_,re_);
volumeout[i] = mag_;
phaseout[i]  = pha_;
} // i
} //  phase	

delete[] dim_sz;
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// to make output volume smaller in size (double -> float)
float *volumeoutF = new float[nelem];
for(i = 0; i < nelem; i++){
volumeoutF[i] = (float)(volumeout[i]);     
}
delete[] volumeout; 

float *volumeoutF_i;
if(phase_flag){
volumeoutF_i = new float[nelem];
for(i = 0; i < nelem; i++){
volumeoutF_i[i] = (float)(phaseout[i]);     
}
delete[] phaseout; 
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// create output nifti image
nifti_image *imout;
if(phase_flag){ imout = nifti_convert_nhdr2nim(spec,nameDataOutMag); }
else{           imout = nifti_convert_nhdr2nim(spec,nameDataOut);    }
free(imout->data);
imout->data = reinterpret_cast<void*> (volumeoutF);
nifti_image_write( imout ); 
nifti_image_free( imout );

nifti_image *phout;
if(phase_flag){
phout = nifti_convert_nhdr2nim(spec,nameDataOutPha);
free(phout->data);
phout->data = reinterpret_cast<void*> (volumeoutF_i);
nifti_image_write( phout ); 
nifti_image_free( phout );
} // phase

return 0;
}
