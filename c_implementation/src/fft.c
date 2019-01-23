#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"

int reverse_matrix_128[]={
	  0,  64,  32,  96,  16,  80,  48, 112,   8,  72,  40, 104,  24,  88,  56, 120,   // 16 each line
	  4,  68,  36, 100,  20,  84,  52, 116,  12,  76,  44, 108,  28,  92,  60, 124,
	  2,  66,  34,  98,  18,  82,  50, 114,  10,  74,  42, 106,  26,  90,  58, 122,
	  6,  70,  38, 102,  22,  86,  54, 118,  14,  78,  46, 110,  30,  94,  62, 126,
	  1,  65,  33,  97,  17,  81,  49, 113,   9,  73,  41, 105,  25,  89,  57, 121,
	  5,  69,  37, 101,  21,  85,  53, 117,  13,  77,  45, 109,  29,  93,  61, 125,
	  3,  67,  35,  99,  19,  83,  51, 115,  11,  75,  43, 107,  27,  91,  59, 123,
	  7,  71,  39, 103,  23,  87,  55, 119,  15,  79,  47, 111,  31,  95,  63, 127
};

int reverse_matrix_64[]={
	  0,  32,  16,  48,  8,   40,  24,  56,
	  4,  36,  20,  52,  12,  44,  28,  60,
	  2,  34,  18,  50,  10,  42,  26,  58,
	  6,  38,  22,  54,  14,  46,  30,  62,
	  1,  33,  17,  49,  9,   41,  25,  57,
	  5,  37,  21,  53,  13,  45,  29,  61,
	  3,  35,  19,  51,  11,  43,  27,  59,
    7,  39,  23,  55,  15,  47,  31,  63
 };
int reverse_matrix_32[]={
	  0,  16,  8,   24,   4,  20,  12,  28,  2,  18,  10,  26,
    6,  22,  14,  30,   1,  17,  9,   25,  5,  21,  13,  29,
	  3,  19,  11,  27,   7,  23,  15,  31
};

void conjugate_complex(int N,complexdouble *in,complexdouble *out)
{
    int i;
	  for(i=0;i<N;i++)
	  {
	  	 out[i].imag = -in[i].imag;
	  	 out[i].real = in[i].real;
	  }
}

void c_abs(complexdouble f[],double out[],int n)
{
	  float t; int i = 0;
	  for(i=0;i<n;i++){
	  	  out[i]= f[i].real * f[i].real+ f[i].imag * f[i].imag;
	  }
}

void c_plus(complexdouble a,complexdouble b,complexdouble *c)
{
	  c->real = a.real + b.real;
	  c->imag = a.imag + b.imag;
}

void c_sub(complexdouble a,complexdouble b,complexdouble *c)
{
	  c->real = a.real - b.real;
	  c->imag = a.imag - b.imag;
}

void c_mul(complexdouble a,complexdouble b,complexdouble *c)
{
	  c->real = a.real * b.real - a.imag * b.imag;
	  c->imag = a.real * b.imag + a.imag * b.real;
}

void c_mul_conjugate(complexdouble a,complexdouble b,complexdouble *c)
{
	  c->real = a.real * b.real - a.imag * b.imag*(-1);
	  c->imag = a.real * (-1)*b.imag + a.imag * b.real;
}

void c_div(complexdouble a,complexdouble b,complexdouble *c)
{
	  c->real = (a.real * b.real + a.imag * b.imag)/(b.real * b.real +b.imag * b.imag);
	  c->imag=(a.imag * b.real - a.real * b.imag)/(b.real * b.real +b.imag * b.imag);
}

void Wn_i(int n,int i,complexdouble *Wn,char flag)
{
	  Wn->real = cos(2*PI*i/n);
	  if(flag == 1)
		    Wn->imag = -sin(2*PI*i/n);
	  else if(flag == 0)
		    Wn->imag = sin(2*PI*i/n);
}

void fft(int N,int M,complexdouble *in)
{
	  complexdouble t,wn;
	  int i,j,k,m,n,l,r,la,lb,lc;
	  for(m=1;m<=M;m++)    //FFT
	  {
	  	  la=pow(2.0,m);
	  	  lb=la/2;
	  	  for(l=1;l<=lb;l++){
	  		    r=(l-1)*pow(2.0,M-m);
	  		    for(n=l-1;n<N-1;n=n+la){
	  			      lc=n+lb;
	  			      Wn_i(N,r,&wn,1);
	  			      c_mul(in[lc],wn,&t);
	  			      c_sub(in[n],t,&(in[lc]));
	  			      c_plus(in[n],t,&(in[n]));
			      }
		    }
	  }
}

void ifft(int N,int M,complexdouble *in)
{
    conjugate_complex(N,in,in);
    fft(N,M,in);
    conjugate_complex(N,in,in);
}

void fft_2D(int mLen,int nLen,int M,int N,complexdouble *A_In,int flag)
{
    FILE *fp; int i,j;
    complexdouble * A;
    int len = mLen * nLen;
    int *b = NULL;
	  if(N==7)
		    b = reverse_matrix_128;  // LUT
	  else if(N==6)
		    b = reverse_matrix_64;
	  else if(N==5)
		    b = reverse_matrix_32;

	  if (flag == DXT_INVERSE){
	      complexdouble *p = A_In;
	  	  for(i=0; i<len; i++)
	  	  {
	  		    (p->imag) *= -1;
	  		    p++;
	  	  }
	  }

	  A = (complexdouble *)malloc(sizeof(complexdouble)*nLen);
	  for(i=0; i<mLen; i++){
		    complexdouble *pr1 = A;
		    complexdouble *pr2 = NULL;
		    complexdouble *pr3 = A_In + (i<<N);
		    for(j=0; j<nLen; j++){
		    	  pr2 = pr3 + b[j];
		    	  pr1->real = pr2->real;
		    	  pr1->imag = pr2->imag;
		    	  pr1++;
		    }
		    fft(nLen,N,A);
		    if (flag == DXT_FORWARD){
			      pr1 = A;
			      pr2 = pr3;
			      for(j=0; j<nLen; j++){
			      	  pr2->real= pr1->real;
			      	  pr2->imag= pr1->imag;
			      	  pr1++;
			      	  pr2++;
			      }
	 	    }
		    else{
		    	  pr1 = A;
		    	  pr2 = pr3;
		    	  for(j=0; j<nLen; j++){
		    		    pr2->real = pr1->real/nLen;
		    		    pr2->imag= pr1->imag/nLen;
		    		    pr1++;
		    		    pr2++;
	          }
        }
    }
    free(A);
    A = (complexdouble *)malloc(sizeof(complexdouble)*mLen);

	  if(M==7)
		    b = reverse_matrix_128;
	  else if(M==6)
		    b = reverse_matrix_64;
	  else if(M==5)
		    b = reverse_matrix_32;

	  for(i=0; i<nLen; i++){
		    complexdouble *pr1 = A;
		    complexdouble *pr2 = NULL;
		    complexdouble *pr3 = A_In + i;
		    for(j=0; j<mLen; j++){
		    	  pr2 = pr3 + (b[j]<<N);
		    	  pr1->real = pr2->real;
		    	  pr1->imag = pr2->imag;
		    	  pr1++;
		    }
		    fft(mLen,M,A);
		    if (flag == DXT_FORWARD){
			      pr1 = A;
			      pr2 = pr3;
			      for(j=0; j<mLen; j++)
			      {
				        pr2->real = pr1->real;
				        pr2->imag = pr1->imag;
				        pr1++;
				        pr2 += nLen;
			      }
		    }
		    else{
			      pr1 = A;
			      pr2 = pr3;
			      for(j=0; j<mLen; j++){
			          pr2->real = pr1->real/mLen;
			          pr2->imag=(-1)*( pr1->imag/mLen);
			          pr1++;
			          pr2 += nLen;
			  }
		}
}
    free(A);
}

void computereal(complexdouble * in,double *out,int length)
{
    int i;
	  for(i=0;i<length;i++)
	  {
		    out[i]=in[i].real;
	  }
}
