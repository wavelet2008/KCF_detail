#ifndef _FHOG_H_
  #define _FHOG_H_
#include <stdio.h>

#include "opencv2/imgproc/imgproc_c.h"
#include "float.h"
// DataType: STRUCT featureMap
// FEATURE MAP DESCRIPTION
// Rectangular map (sizeX x sizeY),
// every cell stores feature vector (dimension = numFeatures)
// map   - matrix of feature vectors to set and get feature vectors (i,j) used formula map[(j * sizeX + i) * p + k], where k - component of feature vector in cell (i, j)
typedef struct{
    int sizeX;
    int sizeY;
    int numFeatures;
    float *map;
} CvLSVMFeatureMapCaskade;


#define PI    CV_PI
#define EPS 0.000001
#define F_MAX FLT_MAX
#define F_MIN -FLT_MAX
// The number of elements in bin
// The number of sectors in gradient histogram building
#define NUM_SECTOR 9
// The number of levels in image resize procedure
// We need Lambda levels to resize image twice
#define LAMBDA 10
#define SIDE_LENGTH 8
#define VAL_OF_TRUNCATE 0.2f
#define LATENT_SVM_OK 0
#define LATENT_SVM_MEM_NULL 2
#define DISTANCE_TRANSFORM_OK 1
#define DISTANCE_TRANSFORM_GET_INTERSECTION_ERROR -1
#define DISTANCE_TRANSFORM_ERROR -2
#define DISTANCE_TRANSFORM_EQUAL_POINTS -3
#define LATENT_SVM_GET_FEATURE_PYRAMID_FAILED -4
#define LATENT_SVM_SEARCH_OBJECT_FAILED -5
#define LATENT_SVM_FAILED_SUPERPOSITION -6
#define FILTER_OUT_OF_BOUNDARIES -7
#define LATENT_SVM_TBB_SCHEDULE_CREATION_FAILED -8
#define LATENT_SVM_TBB_NUMTHREADS_NOT_CORRECT -9
#define FFT_OK 2
#define FFT_ERROR -10
#define LSVM_PARSER_FILE_NOT_FOUND -11

int getFeatureMaps(const IplImage * image, const int k, CvLSVMFeatureMapCaskade **map);
int normalizeAndTruncate(CvLSVMFeatureMapCaskade *map, const float alfa);
int PCAFeatureMaps(CvLSVMFeatureMapCaskade *map);
int allocFeatureMapObject(CvLSVMFeatureMapCaskade **obj, const int sizeX, const int sizeY, const int p);
int freeFeatureMapObject (CvLSVMFeatureMapCaskade **obj);

#endif
