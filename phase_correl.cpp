#include <cmath>
#include <iostream>
#include <cstring>
#include <climits>

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
    void DrawRectangle(unsigned x, unsigned y, unsigned width, unsigned height, unsigned char fill) {
        for (unsigned j = y; j < y + height; j++) {
            for (unsigned i = x; i < x + width; i++) {
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

class PhaseCorrelation {
public:
    PhaseCorrelation() = delete;
    PhaseCorrelation(const PhaseCorrelation&) = delete;
    PhaseCorrelation(PhaseCorrelation&&) = delete;
    PhaseCorrelation &operator=(const PhaseCorrelation&) = delete;
    static void ComputeShift(const GrayscaleImage &image1, const GrayscaleImage &image2, int &deltax, int &deltay) {
        if ((image1.GetWidth() != image2.GetWidth()) || (image1.GetHeight() != image2.GetHeight())) {
            throw std::runtime_error("Wrong image size");
        }

        unsigned width = image1.GetWidth(), height = image1.GetHeight();

        // Convert image pixels to complex number format, use only real part
        double *data1 = new double[(width * height) << 1];
        double *data2 = new double[(width * height) << 1];

        for (unsigned i = 0; i < width * height; i++) {
            data1[i << 1] = (double)image1.GetData()[i];
            data2[i << 1] = (double)image2.GetData()[i];

            data1[(i << 1) + 1] = 0.0;
            data2[(i << 1) + 1] = 0.0;
        } 

        // Perform 2D FFT on each image
        FFT2D(data1, width, height);
        FFT2D(data2, width, height);

        // Compute normalized cross power spectrum
        for (unsigned i = 0; i < width * height; i++) {
            ComputeNormalized(&data1[i << 1], &data2[i << 1], &data1[i << 1]);
        }

        // Perform inversed 2D FFT on obtained matrix
        IFFT2D(data1, width, height);

        for (unsigned i = 0; i < width * height; i++) {
            data1[i << 1] = sqrt(pow(data1[i << 1], 2) + pow(data1[(i << 1) + 1], 2));
        }

        // Search for peak
        double max = 0.0; deltax = 0; deltay = 0;
        for (unsigned j = 0; j < height; j++)
            for (unsigned i = 0; i < width; i++) {
                unsigned offset = i + j * width;
                offset = offset << 1;
                if (data1[offset] > max) {
                    max = data1[offset]; deltax = i; deltay = j;
                }
            }

        delete[] data1;
        delete[] data2;

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
    static void Radix2FFT(double *input, double *output, std::size_t stride) {
        output[0] = input[0] + input[stride << 1];
        output[1] = input[1] + input[(stride << 1) + 1];

        output[2] = input[0] - input[(stride << 1)];
        output[3] = input[1] - input[(stride << 1) + 1];
    }
    static void Sum2FFT(double *input, double *output, std::size_t size, bool inverse) {
        double dfi = (inverse ? 2.0 : -2.0) * M_PI / (double)(size << 1),
            kfi = 0.0;

        for (std::size_t k = 0; k < size; k++) {
            double cosfi = cos(kfi), sinfi = sin(kfi),
                temp[] = {
                    input[k << 1],
                    input[(k << 1) + 1],
                    input[(k + size) << 1],
                    input[((k + size) << 1) + 1]
                };

            output[k << 1] = temp[0] + cosfi * temp[2] - sinfi * temp[3];
            output[(k << 1) + 1] = temp[1] + sinfi * temp[2] + cosfi * temp[3];
            output[(k + size) << 1] = temp[0] - temp[2] * cosfi + temp[3] * sinfi;
            output[((k + size) << 1) + 1] = temp[1] - temp[2] * sinfi - temp[3] * cosfi;
            kfi += dfi;
        }
    }
    static std::size_t *Gen2FFTOffsets(std::size_t size) {
        std::size_t *offsets = new std::size_t[size],
            stride = 2, step;

        std::memset(offsets, 0xff, sizeof(std::size_t) * size);
        offsets[0] = 0;
        if (size > 1) {
            offsets[size >> 1] = 1;
        }

        for (step = size >> 2; step > 0; step >>= 1) {
            std::size_t base = 0;
            for (std::size_t i = 0; i < size; i += step) {
                if (offsets[i] != SIZE_MAX) {
                    base = offsets[i];
                } else {
                    offsets[i] = base + stride;
                }
            }
            stride <<= 1;
        }

        return offsets;
    }
    static void Dit2FFT(double *input, double *output, std::size_t size, bool inverse) {
        std::size_t *offsets = Gen2FFTOffsets(size >> 1),
            stride = 1;

        while (size > 2) {
            stride <<= 1;
            size >>= 1;
        }

        for (std::size_t i = 0; i < stride; i++) {
            Radix2FFT(&input[offsets[i] << 1], &output[i << 2], stride);
        }

        delete[] offsets;

        while (stride > 1) {
            stride >>= 1;
            size <<= 1;
            for (std::size_t i = 0; i < stride; i++) {
                Sum2FFT(&output[(i * size) << 1], &output[(i * size) << 1], size >> 1, inverse);
            }
        }
    }
    static inline void FFT(double *input, double *output, std::size_t size) {
        Dit2FFT(input, output, size, false);
    }
    static inline void IFFT(double *input, double *output, std::size_t size) {
        Dit2FFT(input, output, size, true);
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
                fft_temp[(i + j * width) << 1] = fft_output[j << 1] / static_cast<double>(height);
                fft_temp[((i + j * width) << 1) + 1] = fft_output[(j << 1) + 1] / static_cast<double>(height);
            }
        }
        for (unsigned j = 0; j < height; j++) {
            for (unsigned i = 0; i < width; i++) {
                fft_input[i << 1] = fft_temp[(i + j * width) << 1];
                fft_input[(i << 1) + 1] = fft_temp[((i + j * width) << 1) + 1];
            }
            IFFT(fft_input, fft_output, width);
            for (unsigned i = 0; i < width; i++) {
                input[(i + j * width) << 1] = fft_output[i << 1] / static_cast<double>(width);
                input[((i + j * width) << 1) + 1] = fft_output[(i << 1) + 1] / static_cast<double>(width);
            }
        }

        delete[] fft_input;
        delete[] fft_output;
        delete[] fft_temp;
    }
};

int main()
{
    int deltax, deltay;

    // Generate pair of images
    GrayscaleImage image1(256, 128, 0x00), image2(256, 128, 0xff);
    image1.DrawRectangle(16, 32, 60, 60, 0x80);
    image2.DrawRectangle(8, 40, 60, 60, 0x10);

    PhaseCorrelation::ComputeShift(image1, image2, deltax, deltay);

    std::cout << "Calculated shift: [" << deltax << ", " << deltay << "]" << std::endl;

    return 0;
}