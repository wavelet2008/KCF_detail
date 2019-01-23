#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <time.h>
#include "para.h"
#include "Tracker.h"
#include "fft.h"

int main()
{
    int i,j,h,l,loc,resize_image,scale,heightscale,widthscale,startframe,endframe,vert_delta, horiz_delta,window_sz[2],sz[2],**pos;
	  double *cos_window,*han0,*han1,*fea,*xf,*response,*tmp1,*tmp2,output_sigma,*labels;
	  complexdouble *complex_labels,*complex_xf,*complex_kf,*alphaf,*model_alphaf,*model_xf,*complex_fea,*complex_zf,*complex_tmp1,*complex_tmp2,*complex_kzf;
	  float *imgreyfloat,*imgreyresizefloat;
	  float t1,t2,t3,t4,max;
	  uint8 * im,*rdata,*gdata,*bdata,*imgrey,*imgreyresize,*patch;
	  FILE *fp;
	  clock_t start,end,start1,end1;
	  char img_name[1024]="C:\\Users\\zhang\\KCF_detail\\c_implementation\\data\\ball\\%08d.txt"; char img_path[1024];
	  //init tracking parameters
	  tracker_para* trackerpara =(tracker_para*)malloc(sizeof(tracker_para));
	  target_para*  targetpara  =(target_para*) malloc(sizeof(target_para));
	  para_init(trackerpara);
	  target_para_init(targetpara);

	  im         = (uint8 *)malloc((*targetpara).height*(*targetpara).width*(*targetpara).dim*sizeof(uint8));
	  rdata      = (uint8 *)malloc((*targetpara).height*(*targetpara).width*sizeof(uint8));
	  gdata      = (uint8 *)malloc((*targetpara).height*(*targetpara).width*sizeof(uint8));
	  bdata      = (uint8 *)malloc((*targetpara).height*(*targetpara).width*sizeof(uint8));
	  imgrey     = (uint8 *)malloc((*targetpara).height*(*targetpara).width*sizeof(uint8));
	  imgreyfloat= (float *)malloc((*targetpara).height*(*targetpara).width*sizeof(float));

	  //frame number
	  startframe=(*targetpara).startframe;
	  endframe=(*targetpara).endframe;
	  pos=(int**)malloc((endframe-startframe+1)*sizeof(int*));

	  for(i=0;i<(endframe-startframe+1);i++){
	      pos[i]=(int *)malloc(2*sizeof(int));
	  }
	  resize_image=0; scale=1; output_sigma=0;
	  //read buffer
	  sprintf(img_path,img_name,startframe);
	  readdata(img_path,im,(*targetpara).height*(*targetpara).width*(*targetpara).dim);

    while(sqrt((double)(*targetpara).target_sz[0]*(*targetpara).target_sz[1])>= 100){
	  	  scale=scale*2;
	  	  (*targetpara).target_sz[0] = floor((double)(*targetpara).target_sz[0] / 2);
	  	  (*targetpara).target_sz[1] = floor((double)(*targetpara).target_sz[1] / 2);
	  	  (*targetpara).pos[0]       = floor((double)(*targetpara).pos[0]/2);
	  	  (*targetpara).pos[1]       = floor((double)(*targetpara).pos[1]/2);
	  	  resize_image=1;
	  }

	  //search window size
	  window_sz[0] = floor((double)(*targetpara).target_sz[0] * (1 + (*trackerpara).padding));
	  window_sz[1] = floor((double)(*targetpara).target_sz[1] * (1 + (*trackerpara).padding));

	  //rgb2grey
	  if((*targetpara).dim==3){
	  	  rgb2grey(im,im+(*targetpara).height*(*targetpara).width,im+2*(*targetpara).height*(*targetpara).width,(*targetpara).height,(*targetpara).width,imgrey);
	  }
	  //scale img
	  heightscale       = (int)((*targetpara).height/scale);
	  widthscale        = (int)((*targetpara).width/scale);
	  imgreyresize      = (uint8 *)malloc(heightscale*widthscale*sizeof(uint8));
	  imgreyresizefloat = (float *)malloc(heightscale*widthscale*sizeof(float));

	  if(resize_image){
    		for(i=0;i<(*targetpara).height*(*targetpara).width;i++){
    			imgreyfloat[i]=imgrey[i]*1.0f;
    		}
    		resize_im(imgreyfloat, (*targetpara).width, (*targetpara).height, 1, imgreyresizefloat, widthscale,heightscale);
    		for(i=0;i<heightscale*widthscale;i++){
    			  imgreyresize[i]=(uint8)imgreyresizefloat[i];
    		}
	  }
	  else{
	  	  memcpy(imgreyresize,imgrey,heightscale*widthscale*sizeof(uint8));
	  }
	  output_sigma = sqrt((float)(*targetpara).target_sz[0]*(*targetpara).target_sz[1])*(*trackerpara).output_sigma_factor/(*trackerpara).cell_size;
    //after compute HOG,window size resize
	  sz[0]=floor((double)window_sz[0]/ (*trackerpara).cell_size);
	  sz[1]=floor((double)window_sz[1]/ (*trackerpara).cell_size);

	  labels=(double*)malloc(sz[0]*sz[1]*sizeof(double));
	  gaussian_shaped_labels(output_sigma,sz,labels);
    /*
	  fp=fopen("labels.txt","w");
	  for(i=0; i<sz[0]; i++)
	  {
	  	  for(j=0; j<sz[1]; j++)
	  	  {
	  	  	  fprintf(fp,"%f ",labels[i*sz[1]+j]);
	  	  }
	  	  fprintf(fp,"\n");
	  }
	  fclose(fp); */

	  complex_labels=(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  for(i=0;i<sz[0]*sz[1];i++){
	  	  complex_labels[i].real=labels[i];
	  	  complex_labels[i].imag=0;
	  }
	  //store pre-computed cosine window
	  //fft2d
	  int temp1 = sz[0];
    int temp2 = sz[1];
    int temp3 = (int)(log(sz[0])/log(2));
    int temp4 = (int)(log(sz[1])/log(2));
    fft_2D(temp1, temp2, temp3, temp4,complex_labels, 0);
	  //fft_2D(sz[0], sz[1], (int)(log(sz[0])/log(2)), (int)(log(sz[1])/log(2)),complex_labels, 0);
	  //fft_2D(sz[0],sz[1],log((double)sz[0])/log((double)2),log((double)sz[1])/log((double)2),complex_labels,0);
    /*
	  fp=fopen("complex_labels.txt","w");
	  for(i=0; i<sz[0]; i++){
	  	  for(j=0; j<sz[1]; j++){
	  	  	   fprintf(fp,"%f\n%f\n",complex_labels[i*sz[1]+j].real,complex_labels[i*sz[1]+j].imag);
	  	  }
	  }
	  fclose(fp);*/
	  //store pre-computed cosine window
	  han0=(double*)malloc(sz[0]*sizeof(double));
	  hann(sz[0],han0);
	  han1=(double*)malloc(sz[1]*sizeof(double));
	  hann(sz[1],han1);

    cos_window=(double*)malloc(sz[0]*sz[1]*sizeof(double));
    for(i=0;i<sz[0];i++)
		for(j=0;j<sz[1];j++){
			  cos_window[i*sz[1]+j]=han0[i]*han1[j];
		}

	  //obtain a subwindow for training at newly estimated target position
	  patch=(uint8*)malloc(window_sz[0]*window_sz[1]*sizeof(uint8));
	  get_subwindow(imgreyresize, (*targetpara).pos,window_sz,heightscale,widthscale,patch);

	  xf=(double*)malloc(sz[0]*sz[1]*31*sizeof(double));
	  memset(xf,0,sz[0]*sz[1]*31*sizeof(double));
	  get_features(patch,window_sz[0],window_sz[1],(*trackerpara).cell_size, cos_window,xf);
	  /* sprintf(img_path,"xf%d.txt",startframe);
	  fp=fopen(img_path,"w");
	  for(i=0; i<31; i++){
		    for(j=0;j<sz[0];j++){
		    	  for(h=0;h<sz[1];h++){
		    	      fprintf(fp,"%f ",xf[i*sz[0]*sz[1]+j*sz[1]+h]);
		    	  }
		    	  fprintf(fp,"\n");
		    }
	  }
    fclose(fp);
    */
	  complex_xf=(complexdouble*)malloc(sz[0]*sz[1]*31*sizeof(complexdouble));
	  for(i=0;i<sz[0]*sz[1]*31;i++){
	  	  complex_xf[i].real=xf[i];
	  	  complex_xf[i].imag=0;
	  }
	  for(i=0;i<31;i++){
	  	  fft_2D(sz[0],sz[1],log((double)sz[0])/log((double)2),log((double)sz[1])/log((double)2),complex_xf+i*sz[0]*sz[1],0);
	  }
    /*
	  sprintf(img_path,"complex_xf%d.txt",startframe);
	  fp=fopen(img_path,"w");
	  for(l=0;l<31;l++){
	  	  for(i=0; i<sz[0]; i++){
	  		  for(j=0;j<sz[1];j++){
	  			    fprintf(fp,"%f,%f\n",complex_xf[i*sz[1]+j+l*sz[0]*sz[1]].real,complex_xf[i*sz[1]+j+l*sz[0]*sz[1]].imag);
	  		  }
	  	  }
	  }
	  fclose(fp);*/

	  complex_kf=(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  linear_correlation(complex_xf,complex_xf, sz[0],sz[1],complex_kf);
	  /* sprintf(img_path,"complex_kf%d.txt",startframe);
	  fp=fopen(img_path,"w");
	  for(i=0; i<sz[0]; i++){
	  	  for(j=0;j<sz[1];j++){
	  	  		fprintf(fp,"%f\n",complex_kf[i*sz[1]+j].real);
	  	  		fprintf(fp,"%f\n",complex_kf[i*sz[1]+j].imag);
	  	  }
	  }
	  fclose(fp);*/
	  alphaf=(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  model_alphaf=(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  model_xf=(complexdouble*)malloc(sz[0]*sz[1]*31*sizeof(complexdouble));
	  for(i=0;i<sz[0]*sz[1];i++){
	  	  complex_kf[i].real+= (*trackerpara).lambda;
	  	  c_div(complex_labels[i],complex_kf[i],&alphaf[i]);
	  	  //first frame, train with a single image
	  	  model_alphaf[i] = alphaf[i];
	  }
    /*
	  sprintf(img_path,"model_alphaf%d.txt",startframe);
	  fp=fopen(img_path,"w");
	  for(i=0; i<sz[0]; i++){
	  	  for(j=0;j<sz[1];j++){
	  	  		fprintf(fp,"%f\n%f\n",(float)model_alphaf[i*sz[1]+j].real,(float)model_alphaf[i*sz[1]+j].imag);
	  	  }
	  }
	  fclose(fp);*/
	  for(i=0;i<sz[0]*sz[1]*31;i++){
	  	  model_xf[i] = complex_xf[i];
	  }

	  if(resize_image){
	  	  pos[0][0]=(*targetpara).pos[0]*scale;
	  	  pos[0][1]=(*targetpara).pos[1]*scale;
	  }else{
	  	  pos[0][0]=(*targetpara).pos[0];
	  	  pos[0][1]=(*targetpara).pos[1];
	  }
	  fea=(double*)malloc(sz[0]*sz[1]*31*sizeof(double));
	  complex_fea  =(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  complex_tmp1 =(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  complex_tmp2 =(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  complex_zf   =(complexdouble*)malloc(sz[0]*sz[1]*31*sizeof(complexdouble));
	  complex_kzf  =(complexdouble*)malloc(sz[0]*sz[1]*sizeof(complexdouble));
	  start=clock();
    for(i=startframe+1;i<=endframe;i++){
	      printf("%d frame\n",i);
	  	  start1=clock();
	  	  sprintf(img_path,img_name,i);
	  	  readdata(img_path,im,(*targetpara).height*(*targetpara).width*(*targetpara).dim);
        end1=clock();
		    printf("readdata=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    if((*targetpara).dim==3){
		        rgb2grey(im,im+(*targetpara).height*(*targetpara).width,im+(*targetpara).height*(*targetpara).width*2,(*targetpara).height,(*targetpara).width,imgrey);
		    }
		    end1=clock();
		    printf("rgb2gray=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    if(resize_image){
		    		for(j=0;j<(*targetpara).height*(*targetpara).width;j++){
		    			  imgreyfloat[j]=imgrey[j]*1.0f;
		    		}
		    		resize_im(imgreyfloat, (*targetpara).width, (*targetpara).height, 1, imgreyresizefloat, widthscale,heightscale);
		    		for(j=0;j<heightscale*widthscale;j++){
		    			  imgreyresize[j]=(uint8)imgreyresizefloat[j];
		    		}
		    }
		    else{
		    	  memcpy(imgreyresize,imgrey,heightscale*widthscale*sizeof(uint8));
		    }
		    end1=clock();
		    printf("imresize=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);

		    start1=clock();
		    get_subwindow(imgreyresize, (*targetpara).pos, window_sz,heightscale,widthscale,patch);
		    end1=clock();
		    printf("get_subwindow=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    memset(fea,0,sz[0]*sz[1]*31*sizeof(double));
		    get_features(patch,window_sz[0],window_sz[1],(*trackerpara).cell_size, cos_window,fea);
		    end1=clock();
		    printf("get_features=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    for(j=0;j<sz[0]*sz[1]*31;j++){
		    	  complex_zf[j].real=fea[j];
		    	  complex_zf[j].imag=0;
		    }
		    end1=clock();
		    printf("copydata=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    for(j=0;j<31;j++){
		    	  fft_2D(sz[0],sz[1],log((double)sz[0])/log(2.0),log((double)sz[1])/log(2.0),complex_zf+j*sz[0]*sz[1],0);
		    }
		    end1=clock();
		    printf("fft2d=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();

		    linear_correlation(complex_zf, model_xf,sz[0],sz[1],complex_kzf);
		    end1=clock();
		    printf("correlation=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    for(j=0;j<sz[0]*sz[1];j++){
		    	  c_mul(model_alphaf[j],complex_kzf[j],&complex_tmp2[j]);
		    }
		    fft_2D(sz[0],sz[1],log((double)sz[0])/log(2.0),log((double)sz[1])/log(2.0),complex_tmp2,1);
		    end1=clock();
		    printf("ifft_2d=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    response=(double*)malloc(sz[0]*sz[1]*sizeof(double));
		    computereal(complex_tmp2,response,sz[0]*sz[1]);  //equation for fast detection
		    //target location is at the maximum response. we must take into account the fact that, if the target doesn't move, the peak
		    //will appear at the top-left corner, not at the center (this is discussed in the paper). the responses wrap around cyclically.
		    end1=clock();
		    printf("computereal=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    maxfind(response,sz[0],sz[1],&vert_delta,&horiz_delta);
        if(vert_delta>sz[0]/2)  //wrap around to negative half-space of vertical axis
	          vert_delta = vert_delta - sz[0];
        if(horiz_delta>sz[1]/2)  //same for horizontal axis
	          horiz_delta = horiz_delta -sz[1];
		    end1=clock();
		    printf("maxfind=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    (*targetpara).pos[0]=(*targetpara).pos[0]+(*trackerpara).cell_size*(vert_delta);
		    (*targetpara).pos[1]=(*targetpara).pos[1]+(*trackerpara).cell_size*(horiz_delta);
		    if(resize_image){
		        pos[i-(startframe+1)+1][0]=(*targetpara).pos[0]*scale;
		        pos[i-(startframe+1)+1][1]=(*targetpara).pos[1]*scale;
		    }else{
		        pos[i-(startframe+1)+1][0]=(*targetpara).pos[0];
		        pos[i-(startframe+1)+1][1]=(*targetpara).pos[1];
		    }
		    //obtain a subwindow for training at newly estimated target position
		    end1=clock();
		    printf("computepos=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    get_subwindow(imgreyresize, (*targetpara).pos, window_sz,heightscale,widthscale,patch);

		    end1=clock();
		    printf("get_subwindow=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    memset(fea,0,sz[0]*sz[1]*31*sizeof(double));
		    get_features(patch,window_sz[0],window_sz[1],(*trackerpara).cell_size, cos_window,fea);

		    end1=clock();
		    printf("get_features=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    for(j=0;j<sz[0]*sz[1]*31;j++)
		    {
		    	  complex_xf[j].real=fea[j];
		    	  complex_xf[j].imag=0;
		    }
		    end1=clock();
		    printf("copydata=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    for(j=0;j<31;j++)
		    {
		        fft_2D(sz[0],sz[1],log((double)sz[0])/log(2.0),log((double)sz[1])/log(2.0),complex_xf+j*sz[0]*sz[1],0);
		    }
		    end1=clock();
		    printf("fft_2D=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    //Kernel Ridge Regression, calculate alphas (in Fourier domain)
		    linear_correlation(complex_xf, complex_xf, sz[0],sz[1],complex_kf);
		    end1=clock();
		    printf("correlation=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    for(j=0;j<sz[0]*sz[1];j++){
		    	  complex_kf[j].real+= (*trackerpara).lambda;
		    	  c_div(complex_labels[j],complex_kf[j],&alphaf[j]); //equation for fast training
		    	  //subsequent frames, interpolate model
		    	  model_alphaf[j].real=(1-(*trackerpara).interp_factor)*model_alphaf[j].real+(*trackerpara).interp_factor*alphaf[j].real;
		    	  model_alphaf[j].imag=(1-(*trackerpara).interp_factor)*model_alphaf[j].imag+(*trackerpara).interp_factor*alphaf[j].imag;
		    }
		    /* sprintf(img_path,"model_alphaf%d.txt",i);
		    fp=fopen(img_path,"w");
		    for(j=0; j<sz[0]; j++){
		    		for(h=0;h<sz[1];h++){
		    			  fprintf(fp,"%f\n",model_alphaf[j*sz[1]+h].real);
		    			  fprintf(fp,"%f\n",model_alphaf[j*sz[1]+h].imag);
		    		}
		    }
		    fclose(fp);*/
		    end1=clock();
		    printf("updata_alphaf_model=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    start1=clock();
		    for(j=0;j<sz[0]*sz[1]*31;j++){
		    	  model_xf[j].real=(1-(*trackerpara).interp_factor)*model_xf[j].real;
		    	  model_xf[j].imag=(1-(*trackerpara).interp_factor)*model_xf[j].imag;
		    	  complex_xf[j].real=(*trackerpara).interp_factor*complex_xf[j].real;
		    	  complex_xf[j].imag=(*trackerpara).interp_factor*complex_xf[j].imag;
		    	  c_plus(model_xf[j],complex_xf[j],&model_xf[j]);
		    }
		    end1=clock();
		    printf("updata_xf_model=%lf ms\n",(double)(end1-start1)/CLOCKS_PER_SEC*1000);
		    printf("\n");
	  }
	  end=clock();
	  printf("total=%lf ms\n",(double)(end-start)/CLOCKS_PER_SEC*1000);
	  for(i=0;i<(endframe-startframe+1);i++){
	      printf("%d frame pos=%d,%d\n",i+1,pos[i][0],pos[i][1]);
    }

	  free(trackerpara);
	  free(targetpara);
	  free(im);
	  free(rdata);
    free(gdata);
	  free(bdata);
	  free(imgrey);
	  free(imgreyresize);
	  free(patch);
	  free(pos);
	  free(labels);
	  free(cos_window);
	  free(han0);
	  free(han1);
	  free(fea);
	  free(xf);
	  free(response);
	  free(complex_labels);
	  free(complex_xf);
	  free(complex_kf);
	  free(alphaf);
	  free(model_alphaf);
	  free(model_xf);
	  free(complex_fea);
	  free(complex_zf);
	  free(complex_tmp1);
	  free(complex_tmp2);
	  free(complex_kzf);
    for(i=0;i<(endframe-startframe+1);i++){
		    free(pos[i]);
	  }
	  system("pause");
	  return 0;
}
