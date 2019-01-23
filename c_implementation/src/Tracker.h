#ifndef TRACKER_H
  #define TRACKER_H
  #include "fft.h"
  #include "para.h"

  /*int round(double a);*/
  void para_init(tracker_para* initpara);
  void target_para_init(target_para* initpara);
  int readdata(char *fname,uint8 *im,int datasize);
  int writedata(char *fname,uint8 *im,int datasize);
  void rgb2grey(uint8* rdata,uint8* gdata,uint8* bdata,int height,int width, uint8 *imgGrey);
  void imresize(uint8*in,uint8*out,int height,int width,double scale);
  void circshift(double* in,double *out,int height,int width,int x,int y);
  void gaussian_shaped_labels(double sigma,int* sz,double*labels);
  void hann(int length,double *out);
  void get_subwindow(uint8* im,int* pos,int* sz,int height,int width,uint8* patch);
  void get_features(uint8* im,int height,int width,int cell_size,double* cos_window,double* out);
  void gaussian_correlation(complexdouble *xf,complexdouble *yf,int height,int width,double sigma,complexdouble *kf,int *sz);
  void maxfind(double *response,int height,int width,int* x, int* y);
  void process(double *im,int width,int height,int sbin,double* feat);
  void linear_correlation(complexdouble *xf,complexdouble *yf,int height,int width,complexdouble *xyfsum);
  void alphacopy(float *src, float *dst, struct alphainfo *ofs, int n) ;
  void resize1dtran(float *src, int sheight, float *dst, int dheight, int width, int chan);
  void resize_im(float* src, int sh, int sw, int sc, float* dst, int res_dimy, int res_dimx);

#endif
