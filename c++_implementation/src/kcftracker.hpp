#pragma once

#include "tracker.h"

#ifndef _OPENCV_KCFTRACKER_HPP_
#define _OPENCV_KCFTRACKER_HPP_
#endif

class KCFTracker : public Tracker
{
    public:
        // Constructor
        KCFTracker(bool hog = true, bool fixed_window = true, bool multiscale = true, bool lab = true);
        virtual void init(const cv::Rect &roi, cv::Mat image);
        virtual cv::Rect update(cv::Mat image);

        float interp_factor;       // linear interpolation factor for adaptation
        float sigma;               // gaussian kernel bandwidth
        float lambda;              // regularization
        int   cell_size;           // HOG cell size
        int   cell_sizeQ;          // cell size^2, to avoid repeated operations
        float padding;             // extra area surrounding the target
        float output_sigma_factor; // bandwidth of gaussian target
        int   template_size;       // template size
        float scale_step;          // scale step for multi-scale estimation
        float scale_weight;        // to downweight detection scores of other scales for added stability

    protected:
        // Detect object in the current frame.
        cv::Point2f detect(cv::Mat z, cv::Mat x, float &peak_value);
        // train tracker with a single image
        void train(cv::Mat x, float train_interp_factor);
        // Evaluates a Gaussian kernel with bandwidth SIGMA for all relative shifts between input images X and Y, which must both be MxN. They must    also be periodic (ie., pre-processed with a cosine window).
        cv::Mat gaussianCorrelation(cv::Mat x1, cv::Mat x2);
        // Create Gaussian Peak. Function called only in the first frame.
        cv::Mat createGaussianPeak(int sizey, int sizex);
        // Obtain sub-window from image, with replication-padding and extract features
        cv::Mat getFeatures(const cv::Mat & image, bool inithann, float scale_adjust = 1.0f);
        // Initialize Hanning window. Function called only in the first frame.
        void createHanningMats();
        // Calculate sub-pixel peak for one dimension
        float subPixelPeak(float left, float center, float right);

        cv::Mat _alphaf;
        cv::Mat _prob;
        cv::Mat _tmpl;
        cv::Mat _num;
        cv::Mat _den;
        cv::Mat _labCentroids;

    private:
        int      size_patch[3];
        cv::Mat  hann;
        cv::Size _tmpl_sz;
        float    _scale;
        int      _gaussian_size;
        bool     _hogfeatures;
        bool     _labfeatures;
};
