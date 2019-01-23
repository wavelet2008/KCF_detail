#ifndef __FFT_H__
  #define __FFT_H__
  #include "para.h"
  #define DXT_FORWARD  0
  #define DXT_INVERSE  1
  #define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

  void conjugate_complex(int N,complexdouble *in,complexdouble *out);
  void c_plus(complexdouble a,complexdouble b,complexdouble *c);
  void c_mul(complexdouble a,complexdouble b,complexdouble *c) ;
  void c_mul_conjugate(complexdouble a,complexdouble b,complexdouble *c);
  void c_sub(complexdouble a,complexdouble b,complexdouble *c);
  void c_div(complexdouble a,complexdouble b,complexdouble *c);
  void fft(int N,int M,complexdouble *in);
  void ifft(int N,int M,complexdouble *in);
  void c_abs(complexdouble f[],double out[],int n);
  void computereal(complexdouble * in,double *out,int length);
  void fft_2D(int mLen,int nLen,int M,int N,complexdouble *A_In,int flag);

#endif
