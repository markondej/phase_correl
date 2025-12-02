#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdint>
#include <complex>
#include <vector>

class GrayscaleImage {
public:
    GrayscaleImage(unsigned long width, unsigned long height, unsigned char fill) {
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
    void DrawRectangle(unsigned long x, unsigned long y, unsigned long width, unsigned long height, unsigned char fill) {
        for (unsigned long j = y; j < y + height; j++) {
            for (unsigned long i = x; i < x + width; i++) {
                if ((i < this->width) && (j < this->height)) {
                    data[i + j * this->width] = fill;
                }
            }
        }
    }
    unsigned long  GetWidth() const {
        return width;
    }
    unsigned long  GetHeight() const {
        return height;
    }
    const unsigned char *GetData() const {
        return data;
    }
    void Set(unsigned long width, unsigned long height, const unsigned char *data) {
        delete[] data;
        this->data = new unsigned char[width * height];
        memcpy(this->data, data, sizeof(unsigned char) * width * height);
        this->width = width;
        this->height = height;
    }
private:
    unsigned long  width, height;
    unsigned char *data;
};

class PhaseCorrelation {
public:
    PhaseCorrelation() = delete;
    PhaseCorrelation(const PhaseCorrelation&) = delete;
    PhaseCorrelation(PhaseCorrelation&&) = delete;
    PhaseCorrelation &operator=(const PhaseCorrelation&) = delete;
    static void ComputeShift(const GrayscaleImage &image1, const GrayscaleImage &image2, long &deltax, long &deltay) {
        if ((image1.GetWidth() != image2.GetWidth()) || (image1.GetHeight() != image2.GetHeight()) || !image1.GetWidth() || (image1.GetWidth() & (image1.GetWidth() - 1)) || !image1.GetHeight() || (image1.GetHeight() & (image1.GetHeight() - 1)))
            throw std::runtime_error("Image sizes do not match");

        if (!image1.GetWidth() || (image1.GetWidth() & (image1.GetWidth() - 1)) || !image1.GetHeight() || (image1.GetHeight() & (image1.GetHeight() - 1)))
            throw std::runtime_error("Wrong image size");

        unsigned long width = image1.GetWidth(), height = image1.GetHeight();

        std::vector<std::complex<double>> data1, data2;
        data1.resize(width * height);
        data2.resize(width * height);

        // Convert image pixels to complex number format, use only real part
        for (unsigned long i = 0; i < width * height; i++) {
            data1[i] = { static_cast<double>(image1.GetData()[i]), 0.0 };
            data2[i] = { static_cast<double>(image2.GetData()[i]), 0.0 };
        } 

        // Perform 2D FFT on each image
        FFT2D(data1, width);
        FFT2D(data2, width);

        // Compute normalized cross power spectrum
        for (unsigned long i = 0; i < width * height; i++)
            CompNormCrossCorrelation(data1[i], data2[i], data1[i]);

        // Perform inverse 2D FFT on obtained matrix
        FFT2D(data1, width, true);

        // Search for peak
        std::size_t offset = 0;
        double max = 0.0; deltax = 0; deltay = 0;
        for (unsigned long j = 0; j < height; j++)
            for (unsigned long i = 0; i < width; i++) {
                double d = sqrt(pow(data1[offset].real(), 2) + pow(data1[offset].imag(), 2));
                if (d > max) {
                    max = d;
                    deltax = static_cast<long>(i);
                    deltay = static_cast<long>(j);
                }
                offset++;
            }

        if (deltax > static_cast<long>(width >> 1))
            deltax -= static_cast<long>(width);
        if (deltay > static_cast<long>(height >> 1))
            deltay -= static_cast<long>(height);
    }
private:
    static void CompNormCrossCorrelation(const std::complex<double> &input1, const std::complex<double> &input2, std::complex<double> &output) {
        double diff = std::arg(input1 * std::conj(input2));

        output.real(cos(diff));
        output.imag(sin(diff));
    }
    static void Radix2FFT(const std::complex<double> *input, std::complex<double> *output, std::size_t stride) {
        output[0] = input[0] + input[stride];
        output[1] = input[0] - input[stride];
    }
    static void Sum2FFT(const std::complex<double> *input, std::complex<double> *output, std::size_t size, bool inverse) {
        double dfi = (inverse ? 2.0 : -2.0) * M_PI / static_cast<double>(size << 1),
            kfi = 0.0;

        for (std::size_t k = 0; k < size; k++) {
            double cosfi = cos(kfi), sinfi = sin(kfi);
            std::complex<double> temp[] = { input[k], input[k + size] };

            output[k].real(temp[0].real() + cosfi * temp[1].real() - sinfi * temp[1].imag());
            output[k].imag(temp[0].imag() + sinfi * temp[1].real() + cosfi * temp[1].imag());
            output[k + size].real(temp[0].real() - cosfi * temp[1].real() + sinfi * temp[1].imag());
            output[k + size].imag(temp[0].imag() - sinfi * temp[1].real() - cosfi * temp[1].imag());
            kfi += dfi;
        }
    }
    static std::vector<std::size_t> GenInputOrder(std::size_t size) {
        std::size_t stride = 1;
        std::vector<std::size_t> offset;
        offset.resize(size);

        std::memset(offset.data(), 0xff, sizeof(std::size_t) * size);
        offset[0] = 0;

        for (std::size_t step = size >> 1; step > 0; step >>= 1) {
            std::size_t base = 0;
            for (std::size_t i = 0; i < size; i += step)
                if (offset[i] != SIZE_MAX)
                    base = offset[i];
                else
                    offset[i] = base + stride;
            stride <<= 1;
        }

        return offset;
    }
    static void Dit2FFT(const std::vector<std::complex<double>> &input, std::vector<std::complex<double>> &output, bool inverse = false) {
        std::size_t stride = 1, size = input.size();
        auto offset = GenInputOrder(size >> 1);
        output.resize(size);

        while (size > 2) {
            stride <<= 1;
            size >>= 1;
        }

        for (std::size_t i = 0; i < stride; i++)
            Radix2FFT(&input[offset[i]], &output[i * 2], stride);

        while (stride > 1) {
            stride >>= 1;
            size <<= 1;
            for (std::size_t i = 0; i < stride; i++)
                Sum2FFT(&output[i * size], &output[i * size], size >> 1, inverse);
        }

        if (inverse)
            for (std::size_t i = 0; i < size; i++)
                output[i] /= static_cast<double>(size);
    }
    static void FFT2D(std::vector<std::complex<double>> &data, unsigned long width, bool inverse = false) {
        unsigned long height = data.size() / width;

        auto horizontal_fft = [&]() {
            std::vector<std::complex<double>> input, output;
            input.resize(height);
            for (unsigned long i = 0; i < width; i++) {
                std::size_t offset = i;
                for (unsigned long j = 0; j < height; j++) {
                    input[j] = data[offset];
                    offset += width;
                }
                Dit2FFT(input, output, inverse);
                offset = i;
                for (unsigned long j = 0; j < height; j++) {
                    data[offset] = output[j];
                    offset += width;
                }
            }
        };

        auto vertical_fft = [&]() {
            std::vector<std::complex<double>> input, output;
            input.resize(width);
            for (unsigned long j = 0; j < height; j++) {
                std::size_t offset = j * width;
                for (unsigned long i = 0; i < width; i++) {
                    input[i] = data[offset];
                    offset++;
                }
                Dit2FFT(input, output, inverse);
                offset = j * width;
                for (unsigned long i = 0; i < width; i++) {
                    data[offset] = output[i];
                    offset++;
                }
            }
        };

        if (!inverse) {
            vertical_fft();
            horizontal_fft();
        } else {
            horizontal_fft();
            vertical_fft();
        }
    }
};

int main()
{
    long deltax, deltay;

    try {
        // Generate pair of images
        GrayscaleImage image1(256, 128, 0x00), image2(256, 128, 0xff);
        image1.DrawRectangle(16, 32, 60, 60, 0x80);
        image2.DrawRectangle(8, 40, 60, 60, 0x10);

        PhaseCorrelation::ComputeShift(image1, image2, deltax, deltay);
    } catch (...) {
        std::cerr << "Operation failed" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Calculated shift: [" << deltax << ", " << deltay << "]" << std::endl;
    return EXIT_SUCCESS;
}
