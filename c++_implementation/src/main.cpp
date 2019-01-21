
#include <sstream>
#include <algorithm>
#include "kcftracker.hpp"
#include "dirent.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

cv::Rect_<float> getAxisAlignedBB(std::vector<cv::Point2f> polygon);
std::vector<cv::Rect_<float>> getgroundtruth(std::string txt_file);

int main(int argc, char * argv[])
{
	  bool HOG = true;
	  bool FIXEDWINDOW = false;
	  bool MULTISCALE = false;
	  bool SILENT = true;
	  bool LAB = false;
    std::string sequence = "/sequence";
    if (argc >= 2) {
        sequence = std::string("/vot2014/") + argv[1];
    }
    std::string video_base_path = "D:";
    std::string pattern_jpg = video_base_path + sequence + "/color/*.jpg";
    std::string txt_base_path = video_base_path + sequence + "/groundtruth.txt";
    std::vector<cv::String> image_files;
    cv::glob(pattern_jpg, image_files);
    if (image_files.size() == 0) return -1;
    std::vector<cv::Rect_<float>> groundtruth_rect;
    groundtruth_rect = getgroundtruth(txt_base_path);

    KCFTracker tracker(HOG, FIXEDWINDOW, MULTISCALE, LAB);

    cv::Rect_<float> location = groundtruth_rect[0];
    cv::Mat image; std::vector<cv::Rect_<float>> result_rects;
    int64 tic, toc; double time = 0; bool show_visualization = true;

    for (unsigned int frame = 0; frame < image_files.size(); ++frame) {
        image = cv::imread(image_files[frame]);
        tic = cv::getTickCount();
        if (frame == 0){
          tracker.init( location, image );
        }
        else {
            location = tracker.update(image);
        }
        toc = cv::getTickCount() - tic;
        time += toc;
        result_rects.push_back(location);
        if (show_visualization) {
            cv::putText(image, std::to_string(frame + 1), cv::Point(20, 40), 6, 1,cv::Scalar(0, 255, 255), 2);
            cv::rectangle(image, groundtruth_rect[frame], cv::Scalar(0, 255, 0), 2);
            cv::rectangle(image, location, cv::Scalar(0, 128, 255), 2);
            cv::imshow("KCF", image);
            char key = cv::waitKey(10);
            if (key == 27 || key == 'q' || key == 'Q')
                break;
        }
    }
    time = time / double(cv::getTickFrequency());
    double fps = double(result_rects.size()) / time;
    std::cout << "fps:" << fps << std::endl;
    cv::destroyAllWindows();
    return 0;
}

std::vector<cv::Rect_<float>> getgroundtruth(std::string txt_file)
{
    std::vector<cv::Rect_<float>> rects;
    std::ifstream gt;
    gt.open(txt_file.c_str());
    if (!gt.is_open())
        std::cout << "Ground truth file " << txt_file << " can not be read" << std::endl;
    std::string line; float x1, y1, x2, y2, x3, y3, x4, y4;

    while (getline(gt, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::stringstream ss;
        ss.str(line);
        ss >> x1 >> y1 >> x2 >> y2 >> x3 >> y3 >> x4 >> y4;
        std::vector<cv::Point2f>polygon;
        polygon.push_back(cv::Point2f(x1, y1));
        polygon.push_back(cv::Point2f(x2, y2));
        polygon.push_back(cv::Point2f(x3, y3));
        polygon.push_back(cv::Point2f(x4, y4));
        rects.push_back(getAxisAlignedBB(polygon)); //0-index
    }
    gt.close();
    return rects;
}

cv::Rect_<float> getAxisAlignedBB(std::vector<cv::Point2f> polygon)
{
    double cx = double(polygon[0].x + polygon[1].x + polygon[2].x + polygon[3].x) / 4.;
    double cy = double(polygon[0].y + polygon[1].y + polygon[2].y + polygon[3].y) / 4.;
    double x1 = min(min(min(polygon[0].x, polygon[1].x), polygon[2].x), polygon[3].x);
    double x2 = max(max(max(polygon[0].x, polygon[1].x), polygon[2].x), polygon[3].x);
    double y1 = min(min(min(polygon[0].y, polygon[1].y), polygon[2].y), polygon[3].y);
    double y2 = max(max(max(polygon[0].y, polygon[1].y), polygon[2].y), polygon[3].y);
    double A1 = norm(polygon[1] - polygon[2])*norm(polygon[2] - polygon[3]);
    double A2 = (x2 - x1) * (y2 - y1);
    double s = sqrt(A1 / A2);
    double w = s * (x2 - x1) + 1;
    double h = s * (y2 - y1) + 1;
    cv::Rect_<float> rect(cx-1-w/2.0, cy-1-h/2.0, w, h);
    return rect;
}
