#include <cmath>
#include <iostream>
#include <cstring>
#include <climits>
#include <complex>

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
        std::complex<double> *data1 = new std::complex<double>[(width * height) << 1];
        std::complex<double> *data2 = new std::complex<double>[(width * height) << 1];

        for (unsigned i = 0; i < width * height; i++) {
            data1[i] = { static_cast<double>(image1.GetData()[i]), 0.0 };
            data2[i] = { static_cast<double>(image2.GetData()[i]), 0.0 };
        } 

        // Perform 2D FFT on each image
        FFT2D(data1, width, height);
        FFT2D(data2, width, height);

        // Compute normalized cross power spectrum
        for (unsigned i = 0; i < width * height; i++) {
            ComputeNormalized(data1[i], data2[i], data1[i]);
        }

        // Perform inversed 2D FFT on obtained matrix
        IFFT2D(data1, width, height);

        // Search for peak
        unsigned offset = 0;
        double max = 0.0; deltax = 0; deltay = 0;
        for (unsigned j = 0; j < height; j++)
            for (unsigned i = 0; i < width; i++) {
                double d = sqrt(pow(data1[offset].real(), 2) + pow(data1[offset].imag(), 2));
                if (d > max) {
                    max = d; deltax = i; deltay = j;
                }
                offset++;
            }

        delete[] data1;
        delete[] data2;

        if ((unsigned)deltax > width >> 1)
            deltax = deltax - width;
        if ((unsigned)deltay > height >> 1)
            deltay = deltay - height;
    }
private:
    static void ComputeNormalized(const std::complex<double> &input1, const std::complex<double> &input2, std::complex<double> &output) {
        double a1 = (input1.imag() != 0.0) ? ((input1.real() != 0.0) ? atan(input1.imag() / abs(input1.real())) : std::copysign(M_PI, input1.imag()) / 2.0) : 0.0;
        if (input1.real() < 0.0) {
            a1 = ((input1.imag() < 0.0) ? -1.0 : 1.0) * M_PI - a1;
        }

        double a2 = (input2.imag() != 0.0) ? ((input2.real() != 0.0) ? atan(input2.imag() / abs(input2.real())) : std::copysign(M_PI, input2.imag()) / 2.0) : 0.0;
        if (input2.real() < 0.0) {
            a2 = ((input2.imag() < 0.0) ? -1.0 : 1.0) * M_PI - a2;
        }

        output.real(cos(a1 - a2));
        output.imag(sin(a1 - a2));
    }
    static void Radix2FFT(const std::complex<double> *input, std::complex<double> *output, std::size_t stride) {
        output[0].real(input[0].real() + input[stride].real());
        output[0].imag(input[0].imag() + input[stride].imag());

        output[1].real(input[0].real() - input[stride].real());
        output[1].imag(input[0].imag() - input[stride].imag());
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
    static std::size_t *GenInputOrder(std::size_t size) {
        std::size_t *offset = new std::size_t[size],
            stride = 1;

        std::memset(offset, 0xff, sizeof(std::size_t) * size);
        offset[0] = 0;

        for (std::size_t step = size >> 1; step > 0; step >>= 1) {
            std::size_t base = 0;
            for (std::size_t i = 0; i < size; i += step) {
                if (offset[i] != SIZE_MAX) {
                    base = offset[i];
                } else {
                    offset[i] = base + stride;
                }
            }
            stride <<= 1;
        }

        return offset;
    }
    static void Dit2FFT(const std::complex<double> *input, std::complex<double> *output, std::size_t size, bool inverse = false) {
        std::size_t *offset = GenInputOrder(size >> 1),
            stride = 1;

        while (size > 2) {
            stride <<= 1;
            size >>= 1;
        }

        for (std::size_t i = 0; i < stride; i++) {
            Radix2FFT(&input[offset[i]], &output[i * 2], stride);
        }

        delete[] offset;

        while (stride > 1) {
            stride >>= 1;
            size <<= 1;
            for (std::size_t i = 0; i < stride; i++) {
                Sum2FFT(&output[i * size], &output[i * size], size >> 1, inverse);
            }
        }
    }
    static void FFT2D(std::complex<double> *data, unsigned width, unsigned height) {
        std::complex<double>  *fft_input = new std::complex<double>[width > height ? width : height];
        std::complex<double>  *fft_output = new std::complex<double>[width > height ? width : height];

        for (unsigned j = 0; j < height; j++) {
            unsigned offset = j * width;
            for (unsigned i = 0; i < width; i++) {
                fft_input[i] = data[offset];
                offset++;
            }
            Dit2FFT(fft_input, fft_output, width);
            offset = j * width;
            for (unsigned i = 0; i < width; i++) {
                data[offset] = fft_output[i];
                offset++;
            }
        }
        for (unsigned i = 0; i < width; i++) {
            unsigned offset = i;
            for (unsigned j = 0; j < height; j++) {
                fft_input[j] = data[offset];
                offset += width;
            }
            Dit2FFT(fft_input, fft_output, height);
            offset = i;
            for (unsigned j = 0; j < height; j++) {
                data[offset] = fft_output[j];
                offset += width;
            }
        }

        delete[] fft_input;
        delete[] fft_output;
    }
    static void IFFT2D(std::complex<double> *data, unsigned width, unsigned height) {
        std::complex<double> *fft_input = new std::complex<double>[width > height ? width : height];
        std::complex<double> *fft_output = new std::complex<double>[width > height ? width : height];

        for (unsigned i = 0; i < width; i++) {
            unsigned offset = i;
            for (unsigned j = 0; j < height; j++) {
                fft_input[j] = data[offset];
                offset += width;
            }
            Dit2FFT(fft_input, fft_output, height, true);
            offset = i;
            for (unsigned j = 0; j < height; j++) {
                data[offset] = fft_output[j] / static_cast<double>(height);
                offset += width;
            }
        }
        for (unsigned j = 0; j < height; j++) {
            unsigned offset = j * width;
            for (unsigned i = 0; i < width; i++) {
                fft_input[i] = data[offset];
                offset++;
            }
            Dit2FFT(fft_input, fft_output, width, true);
            offset = j * width;
            for (unsigned i = 0; i < width; i++) {
                data[offset] = fft_output[i] / static_cast<double>(width);
                offset++;
            }
        }

        delete[] fft_input;
        delete[] fft_output;
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