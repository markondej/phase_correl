#include <cmath>
#include <iostream>
#include <cstring>
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#ifdef _WIN32
#ifndef _DEBUG
#pragma comment(lib, "opencv_world3412.lib")
#else
#pragma comment(lib, "opencv_world3412d.lib")
#endif
#endif
#endif

#ifndef USE_OPENCV
class GrayscaleImage {
public:
    GrayscaleImage(unsigned width, unsigned height, unsigned char fill) {
        data = new unsigned char[width * height];
        memset(data, fill, sizeof(unsigned char) * width * height);
        this->width = width;
        this->height = height;
    }
    GrayscaleImage(const GrayscaleImage &source) {
        data = new unsigned char[source.width * source.height];
        memcpy(data, source.data, sizeof(unsigned char) * source.width * source.height);
        width = source.width;
        height = source.height;
    }
    GrayscaleImage(GrayscaleImage&&) = delete;
    GrayscaleImage &operator=(const GrayscaleImage &source) {
        delete[] data;
        data = new unsigned char[source.width * source.height];
        memcpy(data, source.data, sizeof(unsigned char) * source.width * source.height);
        width = source.width;
        height = source.height;
        return *this;
    }
    void DrawRectangle(unsigned startx, unsigned starty, unsigned width, unsigned height, unsigned char fill) {
        for (unsigned j = starty; j < starty + height; j++) {
            for (unsigned i = startx; i < startx + width; i++) {
                if ((i < this->width) && (j < this->height)) {
                    data[i + j * this->width] = fill;
                }
            }
        }
    }
    inline unsigned GetWidth() const {
        return width;
    }
    inline unsigned GetHeight() const {
        return height;
    }
    const unsigned char *GetData() const {
        return data;
    }
    void Set(unsigned width, unsigned height, const unsigned char *data) {
        delete[] data;
        this->data = new unsigned char[width * height];
        memcpy(this->data, data, sizeof(unsigned char) * width * height);
        this->width = width;
        this->height = height;
    }
private:
    unsigned width, height;
    unsigned char *data;
};

#endif
class PhaseCorrelation {
public:
    PhaseCorrelation() = delete;
    PhaseCorrelation(const PhaseCorrelation&) = delete;
    PhaseCorrelation(PhaseCorrelation&&) = delete;
    PhaseCorrelation &operator=(const PhaseCorrelation&) = delete;
#ifndef USE_OPENCV
    static void ComputeShift(const GrayscaleImage &image1, const GrayscaleImage &image2, int &deltax, int &deltay) {
        if ((image1.GetWidth() != image2.GetWidth()) || (image1.GetHeight() != image2.GetHeight())) {
            throw std::runtime_error("Wrong image size");
        }

        unsigned width = image1.GetWidth(), height = image1.GetHeight();
#else
    static void ComputeShift(const cv::Mat &image1, const cv::Mat &image2, int &deltax, int &deltay) {
        if (image1.size != image2.size) {
            throw std::runtime_error("Wrong image size");
        }

        unsigned width = (unsigned)image1.cols, height = (unsigned)image1.rows;
#endif

        // Convert image pixels to complex number format, use only real part
#ifndef USE_OPENCV
        double *data1 = new double[(width * height) << 1];
        double *data2 = new double[(width * height) << 1];

        for (unsigned i = 0; i < width * height; i++) {
            data1[i << 1] = (double)image1.GetData()[i];
            data2[i << 1] = (double)image2.GetData()[i];

            data1[(i << 1) + 1] = 0.0;
            data2[(i << 1) + 1] = 0.0;
        } 
#else
        cv::Mat planes1[] = { cv::Mat_<double>(image1), cv::Mat::zeros(2, image1.size, CV_64F) };
        cv::Mat planes2[] = { cv::Mat_<double>(image2), cv::Mat::zeros(2, image2.size, CV_64F) };
        cv::Mat complex1, complex2;
        cv::merge(planes1, 2, complex1);
        cv::merge(planes2, 2, complex2);
#endif

        // Perform 2D FFT on each image
#ifndef USE_OPENCV
        FFT2D(data1, width, height);
        FFT2D(data2, width, height);
#else
        cv::dft(complex1, complex1);
        cv::dft(complex2, complex2);

        double *data1 = reinterpret_cast<double *>(complex1.data);
        double *data2 = reinterpret_cast<double *>(complex2.data);
#endif

        // Compute normalized cross power spectrum
        for (unsigned i = 0; i < width * height; i++) {
            ComputeNormalized(&data1[i << 1], &data2[i << 1], &data1[i << 1]);
        }

        // Perform inversed 2D FFT on obtained matrix
#ifndef USE_OPENCV
        IFFT2D(data1, width, height);

        for (unsigned i = 0; i < width * height; i++) {
            data1[i << 1] = sqrt(pow(data1[i << 1], 2) + pow(data1[(i << 1) + 1], 2));
        }
#else
        cv::dft(complex1, complex1, cv::DFT_INVERSE);

        cv::split(complex1, planes1);
        cv::magnitude(planes1[0], planes1[1], planes1[0]);
        data1 = reinterpret_cast<double *>(&planes1[0].data[0]);
#endif

        // Search for peak
        double max = 0.0; deltax = 0; deltay = 0;
        for (unsigned j = 0; j < height; j++)
            for (unsigned i = 0; i < width; i++) {
                unsigned offset = i + j * width;
#ifndef USE_OPENCV
                offset = offset << 1;
#endif
                if (data1[offset] > max) {
                    max = data1[offset]; deltax = i; deltay = j;
                }
            }

#ifndef USE_OPENCV
        delete[] data1;
        delete[] data2;

#endif
        if ((unsigned)deltax > width >> 1)
            deltax = deltax - width;
        if ((unsigned)deltay > height >> 1)
            deltay = deltay - height;
    }
private:
    static void ComputeNormalized(double *input1, double *input2, double *output) {
        double a1 = (input1[1] != 0.0f) ? ((input1[0] != 0.0f) ? atan(input1[1] / abs(input1[0])) : M_PI * input1[1] / (2.0f * abs(input1[1]))) : 0.0f;
        if (input1[0] < 0.0f) {
            a1 = ((input1[1] < 0.0f) ? -1.0f : 1.0f) * M_PI - a1;
        }

        double a2 = (input2[1] != 0.0f) ? ((input2[0] != 0.0f) ? atan(input2[1] / abs(input2[0])) : M_PI * input2[1] / (2.0f * abs(input2[1]))) : 0.0f;
        if (input2[0] < 0.0f) {
            a2 = ((input2[1] < 0.0f) ? -1.0f : 1.0f) * M_PI - a2;
        }

        output[0] = cos(a1 - a2);
        output[1] = sin(a1 - a2);
    }
#ifndef USE_OPENCV
    static void DitFFT2(double *input, double *output, unsigned size, int stride, bool invert) {
        if (size > 1) {
            DitFFT2(&input[0], &output[0], size >> 1, stride << 1, invert);
            DitFFT2(&input[stride << 1], &output[size], size >> 1, stride << 1, invert);

            double fi = (invert ? 2.0 : -2.0) * M_PI / (double)size, kfi = 0.0;
            for (unsigned k = 0; k < size; k += 2) {
                double cosfi = cos(kfi), sinfi = sin(kfi), temp[] = {
                    output[k], 
                    output[k + 1], 
                    output[size + k], 
                    output[size + k + 1]
                };

                output[k] = temp[0] + temp[2] * cosfi - temp[3] * sinfi;
                output[k + 1] = temp[1] + temp[2] * sinfi + temp[3] * cosfi;
                output[size + k] = temp[0] - temp[2] * cosfi + temp[3] * sinfi;
                output[size + k + 1] = temp[1] - temp[2] * sinfi - temp[3] * cosfi;
                kfi += fi;
            }
        } else {
            output[0] = input[0];
            output[1] = input[1];
        }
    }
    static inline void FFT(double *input, double *output, unsigned size) {
        DitFFT2(input, output, size, 1, false);
    }
    static inline void IFFT(double *input, double *output, unsigned size) {
        DitFFT2(input, output, size, 1, true);
    }
    static void FFT2D(double *input, unsigned width, unsigned height) {
        double *fft_input = new double[(width > height ? width : height) << 1];
        double *fft_output = new double[(width > height ? width : height) << 1];
        double *fft_temp = new double[(width * height) << 1];

        for (unsigned j = 0; j < height; j++) {
            for (unsigned i = 0; i < width; i++) {
                fft_input[i << 1] = input[(i + j * width) << 1];
                fft_input[(i << 1) + 1] = input[((i + j * width) << 1) + 1];
            }
            FFT(fft_input, fft_output, width);
            for (unsigned i = 0; i < width; i++) {
                fft_temp[(i + j * width) << 1] = fft_output[i << 1];
                fft_temp[((i + j * width) << 1) + 1] = fft_output[(i << 1) + 1];
            }
        }
        for (unsigned i = 0; i < width; i++) {
            for (unsigned j = 0; j < height; j++) {
                fft_input[j << 1] = fft_temp[(i + j * width) << 1];
                fft_input[(j << 1) + 1] = fft_temp[((i + j * width) << 1) + 1];
            }
            FFT(fft_input, fft_output, height);
            for (unsigned j = 0; j < height; j++) {
                input[(i + j * width) << 1] = fft_output[j << 1];
                input[((i + j * width) << 1) + 1] = fft_output[(j << 1) + 1];
            }
        }

        delete[] fft_input;
        delete[] fft_output;
        delete[] fft_temp;
    }
    static void IFFT2D(double *input, unsigned width, unsigned height) {
        double *fft_input = new double[(width > height ? width : height) << 1];
        double *fft_output = new double[(width > height ? width : height) << 1];
        double *fft_temp = new double[(width * height) << 1];

        for (unsigned i = 0; i < width; i++) {
            for (unsigned j = 0; j < height; j++) {
                fft_input[j << 1] = input[(i + j * width) << 1];
                fft_input[(j << 1) + 1] = input[((i + j * width) << 1) + 1];
            }
            IFFT(fft_input, fft_output, height);
            for (unsigned j = 0; j < height; j++) {
                fft_temp[(i + j * width) << 1] = fft_output[j << 1] / (double)height;
                fft_temp[((i + j * width) << 1) + 1] = fft_output[(j << 1) + 1] / (double)height;
            }
        }
        for (unsigned j = 0; j < height; j++) {
            for (unsigned i = 0; i < width; i++) {
                fft_input[i << 1] = fft_temp[(i + j * width) << 1];
                fft_input[(i << 1) + 1] = fft_temp[((i + j * width) << 1) + 1];
            }
            IFFT(fft_input, fft_output, width);
            for (unsigned i = 0; i < width; i++) {
                input[(i + j * width) << 1] = fft_output[i << 1] / (double)width;
                input[((i + j * width) << 1) + 1] = fft_output[(i << 1) + 1] / (double)width;
            }
        }

        delete[] fft_input;
        delete[] fft_output;
        delete[] fft_temp;
    }
#endif
};

int main()
{
    int deltax, deltay;

    // Generate pair of images
#ifndef USE_OPENCV
    GrayscaleImage image1(256, 128, 0x00), image2(256, 128, 0xff);
    image1.DrawRectangle(16, 32, 60, 60, 0x80);
    image2.DrawRectangle(8, 40, 60, 60, 0x10);
#else
    cv::Mat image1(128, 256, CV_8UC1, cv::Scalar(0x00)), image2(128, 256, CV_8UC1, cv::Scalar(0xff));
    cv::rectangle(image1, cv::Rect(16, 32, 60, 60), cv::Scalar(0x80), -1);
    cv::rectangle(image2, cv::Rect(8, 40, 60, 60), cv::Scalar(0x10), -1);
#endif

    PhaseCorrelation::ComputeShift(image1, image2, deltax, deltay);

    std::cout << "Calculated shift: [" << deltax << ", " << deltay << "]" << std::endl;

    return 0;
}