#ifndef PARA_H
  #define MAX(a,b)((a)>(b)?(a):(b))
  #define MIN(a,b)((a)>(b)?(b):(a))

  #define PARA_H
  #define eps 0.0001 // small value, used to avoid division by zero
  #define PI 3.1415926535897932384626433832795028841971

  typedef unsigned char uint8;

  typedef struct tracker_para
  {
  		float padding;
  		float lambda;
  		float output_sigma_factor;
  		float interp_factor;
  		float sigma ;
  		int hog_orientations ;
  		int cell_size ;
  }tracker_para;

  typedef struct target_para{
  	  int width;
  	  int height;
  	  int dim;
  	  int init_rect[4];
  	  int target_sz[2];
  	  int pos[2];
  	  int startframe;
  	  int endframe;
  }target_para;

  typedef struct complexdouble
  {
      double real;
      double imag;
  }complexdouble;

  typedef struct alphainfo {
      int si, di;
      float alpha;
  }alphainfo;
  
#endif
