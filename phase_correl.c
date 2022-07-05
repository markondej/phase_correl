#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

void radix2fft(double *input, double *output, unsigned stride) {
    output[0] = input[0] + input[stride << 1];
    output[1] = input[1] + input[(stride << 1) + 1];

    output[2] = input[0] - input[(stride << 1)];
    output[3] = input[1] - input[(stride << 1) + 1];
}

void sum2fft(double *input, double *output, unsigned size, int inverse) {
    unsigned k;

    double temp[4], cosfi, sinfi,
        dfi = (inverse ? 2.0 : -2.0) * M_PI / (double)(size << 1),
        kfi = 0.0;

    for (k = 0; k < size; k++) {
        cosfi = cos(kfi); sinfi = sin(kfi);
        temp[0] = input[k << 1]; temp[1] = input[(k << 1) + 1];
        temp[2] = input[(k + size) << 1]; temp[3] = input[((k + size) << 1) + 1];

        output[k << 1] = temp[0] + cosfi * temp[2] - sinfi * temp[3];
        output[(k << 1) + 1] = temp[1] + sinfi * temp[2] + cosfi * temp[3];
        output[(k + size) << 1] = temp[0] - temp[2] * cosfi + temp[3] * sinfi;
        output[((k + size) << 1) + 1] = temp[1] - temp[2] * sinfi - temp[3] * cosfi;
        kfi += dfi;
    }
}

unsigned *gen2fftorder(unsigned size) {
    unsigned *offset = (unsigned *)malloc(sizeof(unsigned) * size),
        stride = 1, step;

    memset(offset, 0xff, sizeof(unsigned) * size);
    offset[0] = 0;

    for (step = size >> 1; step > 0; step >>= 1) {
        unsigned base = 0, i;
        for (i = 0; i < size; i += step) {
            if (offset[i] != UINT_MAX) {
                base = offset[i];
            } else {
                offset[i] = base + stride;
            }
        }
        stride <<= 1;
    }

    return offset;
}

void dit2fft(double *input, double *output, unsigned size, int inverse) {
    unsigned *offset = gen2fftorder(size >> 1),
        stride = 1, i;

    while (size > 2) {
        stride <<= 1;
        size >>= 1;
    }

    for (i = 0; i < stride; i++) {
        radix2fft(&input[offset[i] << 1], &output[i << 2], stride);
    }

    free(offset);

    while (stride > 1) {
        stride >>= 1;
        size <<= 1;
        for (i = 0; i < stride; i++) {
            sum2fft(&output[(i * size) << 1], &output[(i * size) << 1], size >> 1, inverse);
        }
    }
}

void fft(double *input, double *output, unsigned size) {
    dit2fft(input, output, size, 0);
}

void ifft(double *input, double *output, unsigned size) {
    dit2fft(input, output, size, 1);
}

void fft2D(double *input, unsigned width, unsigned height) {
    unsigned i, j;

    double *fft_input = (double *)malloc((sizeof(double) * width) << 1);
    double *fft_output = (double *)malloc((sizeof(double) * width) << 1);
    double *fft_temp = (double *)malloc((sizeof(double) * width * height) << 1);

    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            fft_input[i << 1] = input[(i + j * width) << 1];
            fft_input[(i << 1) + 1] = input[((i + j * width) << 1) + 1];
        }
        fft(fft_input, fft_output, width);
        for (i = 0; i < width; i++) {
            fft_temp[(i + j * width) << 1] = fft_output[i << 1];
            fft_temp[((i + j * width) << 1) + 1] = fft_output[(i << 1) + 1];
        }
    }
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            fft_input[j << 1] = fft_temp[(i + j * width) << 1];
            fft_input[(j << 1) + 1] = fft_temp[((i + j * width) << 1) + 1];
        }
        fft(fft_input, fft_output, height);
        for (j = 0; j < height; j++) {
            input[(i + j * width) << 1] = fft_output[j << 1];
            input[((i + j * width) << 1) + 1] = fft_output[(j << 1) + 1];
        }
    }

    free(fft_input);
    free(fft_output);
    free(fft_temp);
}

void ifft2D(double *input, int width, int height) {
    unsigned i, j;

    double *fft_input = (double *)malloc((sizeof(double) * width) << 1);
    double *fft_output = (double *)malloc((sizeof(double) * width) << 1);
    double *fft_temp = (double *)malloc((sizeof(double) * width * height) << 1);

    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            fft_input[j << 1] = input[(i + j * width) << 1];
            fft_input[(j << 1) + 1] = input[((i + j * width) << 1) + 1];
        }
        ifft(fft_input, fft_output, height);
        for (j = 0; j < height; j++) {
            fft_temp[(i + j * width) << 1] = fft_output[j << 1] / (double)height;
            fft_temp[((i + j * width) << 1) + 1] = fft_output[(j << 1) + 1] / (double)height;
        }
    }
    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            fft_input[i << 1] = fft_temp[(i + j * width) << 1];
            fft_input[(i << 1) + 1] = fft_temp[((i + j * width) << 1) + 1];
        }
        ifft(fft_input, fft_output, width);
        for (i = 0; i < width; i++) {
            input[(i + j * width) << 1] = fft_output[i << 1] / (double)width;
            input[((i + j * width) << 1) + 1] = fft_output[(i << 1) + 1] / (double)width;
        }
    }

    free(fft_input);
    free(fft_output);
    free(fft_temp);
}

void computeNormalized(double *f, double *g, double *r) {
    double a1 = (f[1] != 0.0f) ? ((f[0] != 0.0f) ? atan(f[1] / abs(f[0])) : M_PI * f[1] / (2.0f * abs(f[1]))) : 0.0f;
    if (f[0] < 0.0f) {
        a1 = ((f[1] < 0.0f) ? -1.0f : 1.0f) * M_PI - a1;
    }

    double a2 = (g[1] != 0.0f) ? ((g[0] != 0.0f) ? atan(g[1] / abs(g[0])) : M_PI * g[1] / (2.0f * abs(g[1]))) : 0.0f;
    if (g[0] < 0.0f) {
        a2 = ((g[1] < 0.0f) ? -1.0f : 1.0f) * M_PI - a2;
    }

    r[0] = cos(a1 - a2);
    r[1] = sin(a1 - a2);
}

void computeShift(unsigned char *image1, unsigned char *image2, unsigned width, unsigned height, int *deltax, int *deltay) {
    unsigned i, j;

    double* fft_input1 = (double *)malloc((sizeof(double) * width * height) << 1);
    double* fft_input2 = (double *)malloc((sizeof(double) * width * height) << 1);
    double* fft_output = (double *)malloc((sizeof(double) * width * height) << 1);

    // Convert image pixels to complex number format, use only real part
    for (i = 0; i < width * height; i++) {
        fft_input1[i << 1] = (double)image1[i];
        fft_input2[i << 1] = (double)image2[i];

        fft_input1[(i << 1) + 1] = 0.0;
        fft_input2[(i << 1) + 1] = 0.0;
    }

    // Perform 2D FFT on each image
    fft2D(fft_input1, width, height);
    fft2D(fft_input2, width, height);

    // Compute normalized cross power spectrum
    for (i = 0; i < width * height; i++) {
        computeNormalized(&fft_input1[i << 1], &fft_input2[i << 1], &fft_output[i << 1]);
    }

    // Perform inversed 2D FFT on obtained matrix
    ifft2D(fft_output, width, height);

    // Search for peak
    double max = 0.0; *deltax = 0; *deltay = 0;
    for (j = 0; j < height; j++)
        for (i = 0; i < width; i++) {
            double d = sqrt(pow(fft_output[(i + j * width) << 1], 2) + pow(fft_output[((i + j * width) << 1) + 1], 2));
            if (d > max) {
                max = d; *deltax = i; *deltay = j;
            }
        }

    if (*deltax > width >> 1)
        *deltax = *deltax - width;
    if (*deltay > height >> 1)
        *deltay = *deltay - height;

    free(fft_input1);
    free(fft_input2);
    free(fft_output);
}

int main()
{
    unsigned char image1[256 * 128];
    unsigned char image2[256 * 128];

    unsigned i, j;
    int deltax, deltay;

    // Generate pair of images
    for (j = 0; j < 128; j++)
        for (i = 0; i < 256; i++) {
            if ((i >= 16) && (i < 76) && (j >= 32) && (j < 92))
                image1[i + j * 256] = 128;
            else
                image1[i + j * 256] = 0;

            if ((i >= 8) && (i < 68) && (j >= 40) && (j < 100))
                image2[i + j * 256] = 16;
            else
                image2[i + j * 256] = 255;
        }

    computeShift(image1, image2, 256, 128, &deltax, &deltay);

    printf("Calculated shift: [%d, %d]\n", deltax, deltay);

    return EXIT_SUCCESS;
}
