#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "para.h"
#include "Tracker.h"

// unit vectors used to compute gradient orientation
double uu[9] = {1.0000,0.9397,0.7660,0.500,0.1736,-0.1736,-0.5000,-0.7660,-0.9397};
double vv[9] = {0.0000,0.3420,0.6428,0.8660,0.9848,0.9848,0.8660,0.6428,0.3420};

/*int round(double a)
{
    double tmp = a - (int)a;
    if(tmp>=0.5)
        return (int)a + 1;
    else
        return (int)a;
}*/

void para_init(tracker_para* initpara)
{
	  (*initpara).padding             = 1;    //extra area surrounding the target
	  (*initpara).lambda              = 1e-4; //regularization
	  (*initpara).output_sigma_factor = 0.1;  //spatial bandwidth (proportional to target)
	  (*initpara).cell_size           = 4;    //hog
    (*initpara).interp_factor       = 0.02;
    (*initpara).sigma               = 0.5;
    (*initpara).hog_orientations    = 9;
}

void target_para_init(target_para* initpara)
{
	  (*initpara).width        =320;
	  (*initpara).height       =240;
	  (*initpara).init_rect[0] =200;   // left_up x
	  (*initpara).init_rect[1] =114;   // left_up y
	  (*initpara).init_rect[2] =45;    // width
	  (*initpara).init_rect[3] =46;    // height
	  (*initpara).startframe   =1;
	  (*initpara).endframe     =10;
	  (*initpara).dim          =3;
	  (*initpara).target_sz[0] =(*initpara).init_rect[3];
	  (*initpara).target_sz[1] =(*initpara).init_rect[2];
	  (*initpara).pos[0]       =(*initpara).init_rect[1]+floor((double)(*initpara).target_sz[0]/2);
	  (*initpara).pos[1]       =(*initpara).init_rect[0]+floor((double)(*initpara).target_sz[1]/2);
}

int readdata(char *fname,uint8 *im,int datasize)
{
	  FILE *fp = fopen(fname,"r");
	  if (fp==NULL){
	  	  printf("can't open file\n");
	  	  return 1;
	  }
    int  i;
    uint8 * initial_pos = im;
    for(i=0;i<datasize;i++){
      fscanf(fp, "%d", (im+i) );
    }
	  fclose(fp);
	  return 0;
}

int writedata(char *fname,uint8 *im,int datasize)
{
	  FILE *fp = fopen(fname,"wb");
	  if (fp==NULL){
	  	  printf("can't open file\n");
	  	  return 1;
	  }
	  fwrite(im,sizeof(uint8),datasize,fp);
	  fclose(fp);
	  return 0;
}

void rgb2grey(uint8* rdata,uint8* gdata,uint8* bdata,int height,int width, uint8 *imgGrey)
{
	  int pos,tmp;
	  for(pos = 0; pos < width*height; pos++){
	  	  tmp = round(0.2989*rdata[pos] + 0.5870*gdata[pos] + 0.1140*bdata[pos]);   //float point
        tmp = (tmp > 255 ? 255 : tmp);
        tmp = (tmp < 0   ? 0   : tmp);
	  	  imgGrey[pos] = tmp;
	  }
}

void imresize(uint8*in,uint8*out,int height,int width,double scale)
{
	  int i,j,r,c,sq,pos;
	  int row,col;
	  row = (int)(height*scale);
	  col = (int)(width*scale);
	  pos = 0;
	  for (i = 0 ; i<row; i++){
	  	  r= MIN((int)(round((double)i/scale)) ,height-1);
	  	  for (j = 0; j<col; j++){
	  	  		c = MIN((int)(round((double)j/scale)) ,width-1);
	  	  		sq = r * width + c;
	  	  		out[pos++] = in[sq];
	  	  }
	  }
}

void circshift(double* in,double *out,int height,int width,int x,int y)
{
    FILE* fp; int i,j;
    double* tmp=(double *)malloc(width*height*sizeof(double));
    if(x < 0){
		    x=-x;
		    for(i=0;i<(height-x)*width;i++)  //up
		    	  tmp[i]=in[x*width+i];
		    for(i=0;i<x*width;i++)
		    	  tmp[i+(height-x)*width]=in[i];
	  }
	  if(y < 0){
		    y=-y;
		    for(i=0;i<height;i++){         //left
			      for(j=0;j<(width-y);j++)
				        out[i*width+j]=tmp[i*width+j+y];
			      for(j=0;j<y;j++)
				        out[i*width+width-y+j]=tmp[i*width+j];
		    }
	  }
	  free(tmp);
}

void gaussian_shaped_labels(double sigma,int* sz,double*labels)
{
	  FILE* fp; int i,j;
	  int tmp1=(int)floor((double)sz[0]/2);
	  int tmp2=(int)floor((double)sz[1]/2);
	  double*in=(double *)malloc(sz[0]*sz[1]*sizeof(double));
	  memset(in,0,sz[0]*sz[1]*sizeof(double));  //evaluate a Gaussian with the peak at the center element
	  for(i=-8;i<=8;i++)
	  	  for(j=-8;j<=8;j++)
	  		    in[(i+tmp1-1)*sz[1]+(j+tmp2-1)]=exp(-0.5/(sigma*sigma)*(i*i+j*j));

	  tmp1=-tmp1+1;
	  tmp2=-tmp2+1;
	  circshift(in,labels,sz[0],sz[1],tmp1,tmp2);
	  free(in);
}

void hann(int length,double *out)
{
	  int i;
	  for(i=0;i<length;i++)
		    out[i]=0.5*(1-cos(2*PI*i/(length-1)));
}

void get_subwindow(uint8* im,int* pos,int* sz,int height,int width,uint8* patch)
{
	  int i,j; FILE *fp;
	  int *xs=(int *)malloc(sz[0]*sizeof(int));
	  int *ys=(int *)malloc(sz[1]*sizeof(int));
    for(i=0;i<sz[0];i++){
	  	  xs[i]=floor((double)pos[0])+i-floor((double)sz[0]/2);
	  	  if(xs[i]<0)
	  		    xs[i]=0;
	  	  else if(xs[i]>(height-1))
	  		    xs[i]=height-1;
	  }
	  for(i=0;i<sz[1];i++){
	  	  ys[i]=floor((double)pos[1])+i-floor((double)sz[1]/2);
	  	  if(ys[i]<0)
	  		    ys[i]=0;
	      else if(ys[i]>(width-1))
	      	  ys[i]=width-1;
	  }
	  /* fp=fopen("xs.txt","w");
	  for(i=0;i<sz[0];i++)
		    fprintf(fp,"%d\n",xs[i]);
	  fclose(fp);
	  fp=fopen("ys.txt","w");
	  for(i=0;i<sz[1];i++)
		    fprintf(fp,"%d\n",ys[i]);
	  fclose(fp); */
	  for(i=0;i<sz[0];i++)
		  for(j=0;j<sz[1];j++)
			    patch[i*sz[1]+j]=im[xs[i]*width+ys[j]];
	free(xs);  free(ys);
}

void get_features(uint8* im,int height,int width,int cell_size,double* cos_window,double* feat)
{
	  int i,j,h,m,n;
	  double *img,*feat1;
	  FILE *fp;
	  m=height/cell_size-2;
	  n=width/cell_size-2;
	  img=(double*)malloc(height*width*sizeof(double));
	  feat1=(double*)malloc(m*n*32*sizeof(double));
	  for(i=0;i<height*width;i++)
		    img[i]=(double)im[i]/255;
	  process(img,width,height,cell_size,feat1);
	  /*  fp=fopen("fea1.txt","w");
	  for(i=0; i<m*n; i++){
	  		fprintf(fp,"%f\n",feat1[i]);
	  }
	  fclose(fp);*/
	  for(h=0;h<31;h++){
		    for(i=0;i<m;i++){
				    for(j=0;j<n;j++){
						    feat[(i+1)*(n+2)+j+1+h*(m+2)*(n+2)]=feat1[i*n+j+h*m*n];
					  }
			  }
	  }
	  //process with cosine window if needed
	  for(h=0;h<31;h++){
		    for(i=0;i<height/cell_size*width/cell_size;i++){
		        feat[h*height/cell_size*width/cell_size+i]=(feat[h*height/cell_size*width/cell_size+i])*cos_window[i];
		    }
    }
	  /* fp=fopen("feat3.txt","w");
	  for(i=0; i<(m+2); i++){
	  	  for(j=0;j<(n+2);j++)
	  		    fprintf(fp,"%f  ",feat[i*(n+2)+j]);
	  	      fprintf(fp,"\n");
	  }
	  fclose(fp); */
	  free(img); free(feat1);
}

void linear_correlation(complexdouble *xf,complexdouble *yf,int height,int width,complexdouble *xyfsum)
{
	  FILE* fp;  int i,j;
	  int N=height*width;
	  complexdouble *xyf=(complexdouble *)malloc(N*31*sizeof(complexdouble));
	  memset(xyfsum,0,N*sizeof(complexdouble));
	  for(i=0;i<N*31;i++){
	  	  c_mul_conjugate(xf[i],yf[i],&xyf[i]);
	  }
	  /* fp=fopen("xyf.txt","w");
	  for(i=0;i<N*31;i++){
	  	fprintf(fp,"%f\n",xyf[i].real);
	  }
	  fclose(fp); */
	  for(i=0;i<31;i++){
	  	  for(j=0;j<N;j++){
	  	  	  xyfsum[j].real+=xyf[i*N+j].real;
	  	  	  xyfsum[j].imag+=xyf[i*N+j].imag;
	  	  }
	  }
	  /* fp=fopen("xyfsum1.txt","w");
	  for(i=0;i<N;i++){
	  	  fprintf(fp,"%f\n",xyfsum[i].real);
	  }
	  fclose(fp); */
	  for(i=0;i<N;i++){
	  		xyfsum[i].real=xyfsum[i].real/(N*31);
	  		xyfsum[i].imag=xyfsum[i].imag/(N*31);
	  }
	  /* fp=fopen("xyfsum2.txt","w");
	  for(i=0;i<N;i++){
	  	  fprintf(fp,"%f\n",xyfsum[i].real);
	  }
	  fclose(fp); */
	  free(xyf);
}

void maxfind(double *response,int height,int width,int* x, int* y)
{
	  int i,j;
	  double max=-1e10;
	  for(i=0;i<height;i++)
	  	  for(j=0;j<width;j++)
	  		{
	  			  if(response[i*width+j]>max)
	  			  {
	  				    max=response[i*width+j];
	  				    *x=i;
	  				    *y=j;
	  			  }
	  		}
}
// main function: takes a double color image and a bin size and returns HOG features
void process(double *im,int width,int height,int sbin,double* feat)
{
	  int dims[2];
	  // memory for caching orientation histograms & their norms
	  int blocks[2];
	  // memory for HOG features
	  int out[3]; int visible[2];
	  int i,j,h,x,y,best_o,o,ixp,iyp;
	  double *s,dx,dy,v,best_dot,dot,xp,yp,vx0,vy0,vx1,vy1,*src1,*src2,*dst,*end,*src, *p, n1, n2, n3, n4,t1,t2,t3,t4,h1,h2,h3,h4,sum,*mag,tmp1,tmp2,tmp3,tmp4;
	  double *hist,*norm; 	FILE *fp;

	  dims[0]=width;
	  dims[1]=height;
	  blocks[0] = (int)round((double)dims[0]/(double)sbin);
	  blocks[1] = (int)round((double)dims[1]/(double)sbin);
	  hist = (double *)malloc(blocks[0]*blocks[1]*18*sizeof(double));
	  norm = (double *)malloc(blocks[0]*blocks[1]*sizeof(double));
	  memset(hist,0,blocks[0]*blocks[1]*18*sizeof(double));

	  out[0]     = MAX(blocks[0]-2, 0);
	  out[1]     = MAX(blocks[1]-2, 0);
	  out[2]     = 27+4+1;
	  visible[0] = blocks[0]*sbin;
	  visible[1] = blocks[1]*sbin;
	  mag=(double *)malloc((visible[0]-2)*(visible[1]-2)*sizeof(double));
    fp=fopen("hog.txt","w");
	  for (x = 1; x < visible[1]-1; x++){
		    for (y = 1; y < visible[0]-1; y++){
		    	  s = im + MIN(x, dims[1]-2)*dims[0] + MIN(y, dims[0]-2);  // first color channel
		    	  dy = *(s+1) - *(s-1);
		    	  dx = *(s+dims[0]) - *(s-dims[0]);
		    	  v = dx*dx + dy*dy; // snap to one of 18 orientations
		    	  best_dot = 0;
		    	  best_o = 0;
		    	  for (o = 0; o < 9; o++){
		    		    dot = uu[o]*dx + vv[o]*dy;
		    		  if (dot > best_dot){
		    		  	  best_dot = dot;
		    		  	  best_o = o;
		    		  }
		    		  else if (-dot > best_dot){
		    		  	best_dot = -dot;
		    		  	best_o   = o + 9;
		    		  }
		    	  }
		    	  // add to 4 histograms around pixel using linear interpolation
		    	  xp  = ((double)x+0.5)/(double)sbin - 0.5;
		  	    yp  = ((double)y+0.5)/(double)sbin - 0.5;
		  	    ixp = (int)floor(xp);
		  	    iyp = (int)floor(yp);

			      vx0 = xp-ixp;
			      vy0 = yp-iyp;
			      vx1 = 1.0-vx0;
			      vy1 = 1.0-vy0;

			      v = sqrt(v);
			      fprintf(fp,"%f\n",v);
			      if (ixp >= 0 && iyp >= 0){
			      	  *(hist + ixp*blocks[0] + iyp + best_o*blocks[0]*blocks[1])+=vx1*vy1*v;
			      }
			      if (ixp+1 < blocks[1] && iyp >= 0){
			      	  *(hist + (ixp+1)*blocks[0] + iyp + best_o*blocks[0]*blocks[1])+=vx0*vy1*v;
			      }
			      if (ixp >= 0 && iyp+1 < blocks[0]){
			      	  *(hist + ixp*blocks[0] + (iyp+1) + best_o*blocks[0]*blocks[1])+=vx1*vy0*v;
			      }
			      if (ixp+1 < blocks[1] && iyp+1 < blocks[0]){
			      	  *(hist + (ixp+1)*blocks[0] + (iyp+1) + best_o*blocks[0]*blocks[1])+=vx0*vy0*v;
			      }
        }
    }
    fclose(fp);
	  fp=fopen("hoghist.txt","w");
    for(i=0;i<blocks[0]*blocks[1]*18;i++)
        fprintf(fp,"%f\n",hist[i]);
    fclose(fp);
    memset(norm,0,blocks[0]*blocks[1]*sizeof(double));
	  // compute energy in each block by summing over orientations
	  for (o = 0; o < 9; o++)
	  {
		    src1 = hist + o*blocks[0]*blocks[1];
		    src2 = hist + (o+9)*blocks[0]*blocks[1];
		    dst = norm;
		    end = norm + blocks[1]*blocks[0];
		    while (dst < end)
		    {
		    	  *(dst++) += (*src1 + *src2) * (*src1 + *src2);
		    	  src1++;
		    	  src2++;
		    }
	  }
	  /* fp=fopen("hognorm.txt","w");
    for(i=0;i<blocks[0]*blocks[1];i++)
        fprintf(fp,"%f\n",norm[i]);
    fclose(fp);
    */
	  fp=fopen("p.txt","w");
	  // compute features
	  for (x = 0; x < out[1]; x++)
	  {
		    for (y = 0; y < out[0]; y++)
		    {
			      dst = feat + x*out[0] + y;
			      p = norm + (x+1)*blocks[0] + y+1;
			      tmp1= *p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps;
			      tmp1= round(tmp1*10e3)/10e3;
			      n1 = 1.0 / sqrt(tmp1);
			      p = norm + (x+1)*blocks[0] + y;
			      tmp2= *p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps;
			      tmp2= round(tmp2*10e3)/10e3;
			      n2 = 1.0 / sqrt(tmp2);
			      p = norm + x*blocks[0] + y+1;
			      tmp3=  *p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps;
			      tmp3= round(tmp3*10e3)/10e3;
			      n3 = 1.0 / sqrt(tmp3);
			      p = norm + x*blocks[0] + y;
			      tmp4=*p + *(p+1) + *(p+blocks[0]) + *(p+blocks[0]+1) + eps;
			      tmp4=round(tmp4*10e3)/10e3;
			      n4 = 1.0 / sqrt(tmp4);
		        fprintf(fp,"%f,%f,%f,%f\n",n1,n2,n3,n4);
			      t1 = 0;
			      t2 = 0;
			      t3 = 0;
			      t4 = 0;
			      // contrast-sensitive features
			      src = hist + (x+1)*blocks[0] + (y+1);
			      for (o = 0; o < 18; o++){
				        h1   = MIN(*src * n1, 0.2);
				        h2   = MIN(*src * n2, 0.2);
				        h3   = MIN(*src * n3, 0.2);
				        h4   = MIN(*src * n4, 0.2);
                *dst = 0.5 * (h1 + h2 + h3 + h4);
				        t1  += h1;
				        t2  += h2;
				        t3  += h3;
				        t4  += h4;
				        dst += out[0]*out[1];
				        src += blocks[0]*blocks[1];
		    	  }
			      src = hist + (x+1)*blocks[0] + (y+1); // contrast-insensitive features
			      for (o = 0; o < 9; o++)
			      {
			      	  sum  = *src + *(src + 9*blocks[0]*blocks[1]);
			      	  h1   = MIN(sum * n1, 0.2);
			      	  h2   = MIN(sum * n2, 0.2);
			      	  h3   = MIN(sum * n3, 0.2);
			      	  h4   = MIN(sum * n4, 0.2);
			      	  *dst = 0.5 * (h1 + h2 + h3 + h4);
			      	  dst += out[0]*out[1];
			      	  src += blocks[0]*blocks[1];
			      }
			      *dst = 0.2357 * t1;  // texture features
			      dst += out[0]*out[1];
			      *dst = 0.2357 * t2;
			      dst += out[0]*out[1];
			      *dst = 0.2357 * t3;
			      dst += out[0]*out[1];
			      *dst = 0.2357 * t4;
			      dst += out[0]*out[1]; // truncation feature
			      *dst = 0;
		    }
	  }
	  fclose(fp); free(mag); free(hist); free(norm);
}
void gaussian_correlation(complexdouble *xf,complexdouble *yf,int height,int width,double sigma,complexdouble *xyfsum,int *sz)
{
    FILE* fp; int i,j;  double sx=0; double sy=0;
    int N=height*width;
    double *xx   =(double *)malloc(N*31*sizeof(double));
    double *yy   =(double *)malloc(N*31*sizeof(double));
    double *xy   =(double *)malloc(N*31*sizeof(double));
    double *xysum=(double *)malloc(N*sizeof(double));
    complexdouble *xyf=(complexdouble *)malloc(N*31*sizeof(complexdouble));
    memset(xysum,0,N*sizeof(double));
    c_abs(xf,xx,N*31);
    c_abs(yf,yy,N*31);
    /* fp=fopen("xx.txt","w");
    for(i=0;i<N*31;i++)
		    fprintf(fp,"%f\n",xx[i]);
	  fclose(fp);
	  fp=fopen("yy.txt","w");
	  for(i=0;i<N*31;i++)
	  	  fprintf(fp,"%f\n",yy[i]);
	  fclose(fp); */
	  for(i=0;i<N*31;i++){
		    sx+=xx[i];
		    sy+=yy[i];
	  }
	  printf("%f,%f\n",sx,sy);
	  sx/=N;
	  sy/=N;
	  printf("%f,%f\n",sx,sy);
	  /* fp=fopen("xf.txt","w");
	  for(i=0;i<N;i++)
		    fprintf(fp,"%f,%f ",xf[i].real,xf[i].imag);
	  fclose(fp);
	  fp=fopen("yf.txt","w");
	  for(i=0;i<N;i++)
		    fprintf(fp,"%f,%f ",yf[i].real,yf[i].imag);
	  fclose(fp); */
	  for(i=0;i<N*31;i++){
	  	  c_mul_conjugate(xf[i],yf[i],&xyf[i]);
	  }
	  /* fp=fopen("xyf.txt","w");
	  for(i=0;i<N;i++)
	  	  fprintf(fp,"%f,%f\n",xyf[i].real,xyf[i].imag);
	  fclose(fp); */
	  for(i=0;i<31;i++){
	  	  fft_2D(height,width,log((double)height)/log(2.0),log((double)width)/log(2.0),xyf+i*height*width,1);
	  }
	  computereal(xyf,xy,N*31);
	  for(i=0;i<31;i++){
	  	  for(j=0;j<N;j++)
	  		    xysum[j]+=xy[i*N+j];
	  }
	  for(i=0;i<N;i++){   //cross-correlation term in Fourier domain
	  	  xyfsum[i].real=exp(-1/(sigma*sigma)*MAX(0,(sx+sy-2*xysum[i])/(N*31)));
	  	  xyfsum[i].imag=0;
	  }
    //calculate gaussian response for all positions, then go back to the Fourier domain
	  fft_2D(height,width,log((double)height)/log(2.0),log((double)width)/log(2.0),xyfsum,0);
	  free(xx); free(yy); free(xyf); free(xy);
}
// copy src into dst using precomputed interpolation values
void alphacopy(float *src, float *dst, struct alphainfo *ofs, int n)
{
    struct alphainfo *end = ofs + n;
    while (ofs != end)
	  {
        dst[ofs->di] += ofs->alpha * src[ofs->si];
        ofs++;
    }
}
// resize along each column, result is transposed, so we can apply it twice for a complete resize
void resize1dtran(float *src, int sheight, float *dst, int dheight, int width, int chan)
{
	  int c,x; int dy,sy; int k = 0;
    float scale = (float)dheight/(float)sheight;
    float invscale = (float)sheight/(float)dheight;
    int len = (int)ceil(dheight*invscale)+2*dheight;  //we cache the interpolation values since they can be shared among different columns
	  alphainfo* ofs=(alphainfo*)malloc(len*sizeof(alphainfo)); //struct alphainfo ofs[len];

    for (dy = 0; dy < dheight; dy++){
        float fsy1 = dy * invscale;
        float fsy2 = fsy1 + invscale;
        int sy1 = (int)ceil(fsy1);
        int sy2 = (int)floor(fsy2);
        if (sy1 - fsy1 > 1e-3){
            ofs[k].di = dy*width;
            ofs[k].si = sy1-1;
            ofs[k++].alpha = (sy1 - fsy1) * scale;
        }
        for (sy = sy1; sy < sy2; sy++){
            ofs[k].di = dy*width;
            ofs[k].si = sy;
            ofs[k++].alpha = scale;
        }
        if (fsy2 - sy2 > 1e-3){
            ofs[k].di = dy*width;
            ofs[k].si = sy2;
            ofs[k++].alpha = (fsy2 - sy2) * scale;
        }
    }
    // resize each column of each color channel
	  memset(dst,0,chan*width*dheight*sizeof(float));
    for (c = 0; c < chan; c++){
        for (x = 0; x < width; x++){
            float *s = src + c*width*sheight + x*sheight;
            float *d = dst + c*width*dheight + x;
            alphacopy(s, d, ofs, k);
        }
    }
	  free(ofs);
}

// The rescaled image has the same channels as the original image
void resize_im(float* src, int sh, int sw, int sc, float* dst, int res_dimy, int res_dimx)
{
	  float* rescaledim = dst;
	  float* tempim = (float*)malloc(sizeof(float) * sc*sw*res_dimy);
	  resize1dtran(src, sh, tempim, res_dimy, sw, sc);
	  resize1dtran(tempim, sw, rescaledim, res_dimx, res_dimy, sc);
	  free(tempim);
}
