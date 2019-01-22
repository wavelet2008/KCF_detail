# KCF_detail（Kenelized Correlation Filter设计思想及TI DSP C6678实现）

## idea（理念）

* 1.原理
  * 论文
  * 知乎上的相关blog
* 2.c++实现版
* 3.python实现版
* 4.matlab实现版
* 5.嵌入式移植
  * TI C6678
  * Xilinx Virtex 7 
  * GPU TX2 ( support later..）

## 技术栈

* c++ & c
* matlab
* python
* FFT
* HOG
* verilog


## 研究方向（在此处贴链接）

* [KCF作者网站](http://www.robots.ox.ac.uk/~joao/)
* [鲤鱼王KCF知乎博客](https://zhuanlan.zhihu.com/p/26685032)
* [Knight知乎博客KCF](https://zhuanlan.zhihu.com/p/33543297)
* [KCF公式推导错误及验证YaqiLYU](https://zhuanlan.zhihu.com/p/26766703)

## 文件路径

* paper_blog

* c++__implementation

* matlab_implementation

* python_implementation

* dsp_implementation


## 移植工作
* opencv
  * opencv  MAT、CvMat、IplImage类实现
  * opencv  HOG实现（c++ 已有算法细节）
  * opencv  FFT实现
* TI library DSPLIB
  ![dsplib  fft](https://zhangxin-1258487717.cos.ap-beijing.myqcloud.com/dsplib_fft.jpg) 
* [emcv一个2012年的dsp opencv库](https://github.com/zhangxin6/emcv)

## C++ implementation 解析

### 文件架构

![kcf_c++](https://zhangxin-1258487717.cos.ap-beijing.myqcloud.com/kcf_detail_file_dir.jpg) 

**kcf c++代码汇总表**
{% raw %}

<table>
  <tr>   
    <td colspan="2">文件名称</td>
    <td>行数</td>
    <td>描述</td>    
    <td>函数</td>
  </tr>
  <tr>
    <td rowspan="3">.cpp</td>
    <td>main.cpp</td>
    <td>115</td>
    <td>主程序</td>
    <td>
       int main(int argc, char * argv[]) 函数入口，主程序，54行代码<br>
       std::vector&lt;cv::Rect_&lt;float&gt;&gt; getgroundtruth(std::string txt_file) 获取目标位置真值，24行<br>
       cv::Rect_&lt;float&gt; getAxisAlignedBB(std::vector&lt;cv::Point2f&gt; polygon) 获取一个坐标轴对齐后的矩形框，15行
    </td>
  </tr>
  <tr>
  <td>kcf_tracker.cpp</td>
  <td>351</td>
  <td>kcf跟踪器</td>
  <td>
     KCFTracker::KCFTracker(bool hog, bool fixed_window, bool multiscale, bool lab) 跟踪器构造函数39行<br>
     void KCFTracker::init(const cv::Rect &amp;roi, cv::Mat image) 跟踪器初始化54行<br>
     cv::Rect KCFTracker::update(cv::Mat image) 跟踪器更新位置44行 <br>
     cv::Point2f KCFTracker::detect(cv::Mat z, cv::Mat x, float &amp;peak_value) 检测21行<br>
     void KCFTracker::train(cv::Mat x, float train_interp_factor) 训练15行<br>
     cv::Mat KCFTracker::gaussianCorrelation(cv::Mat x1, cv::Mat x2) 高斯相关29行<br>
     cv::Mat KCFTracker::createGaussianPeak(int sizey, int sizex) 高斯峰值获取16行<br>
     cv::Mat KCFTracker::getFeatures(const cv::Mat &amp; image, bool inithann, float scale_adjust) 特征获取120行<br>
     void KCFTracker::createHanningMats() 初始化汉明窗22行<br>
     float KCFTracker::subPixelPeak(float left, float center, float right) 一维子峰值计算10行<br>
  </td>
  </tr>
  <tr>
  <td>fhog.cpp</td>
  <td>326</td>
  <td>HOG特征提取</td>
  <td>
     int getFeatureMaps(const IplImage* image, const int k, CvLSVMFeatureMapCaskade **map) 获取特征图126行<br>
     int normalizeAndTruncate(CvLSVMFeatureMapCaskade *map, const float alfa) 归一化及截断，85行<br>
     int PCAFeatureMaps(CvLSVMFeatureMapCaskade *map) 特征降维，52行<br>
     int allocFeatureMapObject(CvLSVMFeatureMapCaskade **obj, const int sizeX,const int sizeY, const int numFeatures) 分配特征表14行<br>
     int freeFeatureMapObject (CvLSVMFeatureMapCaskade **obj) 释放特征表空间8行<br>
  </td>
  </tr>
  <tr>
    <td rowspan="4">.hpp</td>
    <td>kcftracker.hpp</td>
    <td>60</td>
    <td>kcf跟踪器头文件</td>
    <td>
      KCFTracker类，无子函数<br>
    </td>
 </tr>
 <tr>		
    <td>ffttools.hpp</td>
    <td>179</td>
    <td>fft头文件</td>
    <td>
       cv::Mat fftd(cv::Mat img, bool backwards) 快速FFT，52行代码<br>
       cv::Mat real(cv::Mat img) 获取复数实部，6行代码<br>
       cv::Mat imag(cv::Mat img) 获取复数虚部，6行代码<br>
       cv::Mat magnitude(cv::Mat img) 获取复数模值，10行代码<br>
       cv::Mat complexMultiplication(cv::Mat a, cv::Mat b) 复数乘法，15行代码<br>
       cv::Mat complexDivision(cv::Mat a, cv::Mat b) 复数除法，18行代码<br>
       void rearrange(cv::Mat &amp;) 图像重排，19行代码<br>
       void normalizedLogTransform(cv::Mat &amp;img)  归一化log变换，7行代码<br>
    </td>
  </tr>
  <tr>
    <td>rectools.hpp</td>
    <td>91</td>
    <td>rect头文件</td>
    <td>
       inline cv::Vec&lt;t, 2 &gt; center(const cv::Rect_&lt;t&gt; &amp;rect) 获取中心 4行<br>
       inline t x2(const cv::Rect_&lt;t&gt; &amp;rect) 4行<br>
       inline t y2(const cv::Rect_&lt;t&gt; &amp;rect) 4行<br>
       inline void resize(cv::Rect_&lt;t&gt; &amp;rect, float scalex, float scaley = 0) 调整尺寸10行<br>
       inline void limit(cv::Rect_&lt;t&gt; &amp;rect, cv::Rect_&lt;t&gt; limit) 边界处理17行<br> 
       inline void limit(cv::Rect_&lt;t&gt; &amp;rect, t width, t height, t x = 0, t y = 0) 3行<br>
       inline cv::Rect getBorder(const cv::Rect_&lt;t &gt; &amp;original, cv::Rect_&lt;t &gt; &amp; limited) 获取边界10行<br>
       inline cv::Mat subwindow(const cv::Mat &amp;in, const cv::Rect &amp; window, int borderType = cv::BORDER_CONSTANT) 子窗口14行<br>
       inline cv::Mat getGrayImage(cv::Mat img) 获取灰度图6行<br> 
    </td>
  </tr>
  <tr>
     <td>lambda.hpp</td>
    <td>18</td>
    <td>一些预定义参数</td>
    <td>
       无子函数<br>
    </td>
  </tr>
  <tr>
    <td rowspan="2">.h</td>
    <td>tracker.h</td>
    <td>16</td>
    <td>tracker class 定义</td>
    <td>
       无子函数<br>
    </td>
  </tr> 
  <tr>
    <td>dirent.h</td>
    <td>1225</td>
    <td>linux 底层在windows上的实现</td>
    <td>
       过长，不展开<br>
    </td>
  </tr>  
</table>
{% endraw %}

