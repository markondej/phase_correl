#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void ditfft2(double *input, double *output, unsigned size, unsigned stride, int invert) {
    double temp[4], fi, kfi, cosfi, sinfi;
    unsigned k;

    if (size > 1) {
        ditfft2(&input[0], &output[0], size >> 1, stride << 1, invert);
        ditfft2(&input[stride << 1], &output[size], size >> 1, stride << 1, invert);

        if (!invert) {
            fi = -2.0 * M_PI / (double)size; kfi = 0.0;
        } else {
            fi = 2.0 * M_PI / (double)size; kfi = 0.0;
        }
        for (k = 0; k < size; k += 2) {
            cosfi = cos(kfi); sinfi = sin(kfi);
            temp[0] = output[k]; temp[1] = output[k + 1];
            temp[2] = output[size + k]; temp[3] = output[size + k + 1];

            output[k] = temp[0] + temp[2] * cosfi - temp[3] * sinfi;
            output[k + 1] = temp[1] + temp[2] * sinfi + temp[3] * cosfi;
            output[size + k] = temp[0] - temp[2] * cosfi + temp[3] * sinfi;
            output[size + k + 1] = temp[1] - temp[2] * sinfi - temp[3] * cosfi;
            kfi += fi;
        }
    }
    else {
        output[0] = input[0];
        output[1] = input[1];
    }
}

void fft(double *input, double *output, unsigned size) {
    ditfft2(input, output, size, 1, 0);
}

void ifft(double *input, double *output, unsigned size) {
    ditfft2(input, output, size, 1, 1);
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
    double a1 = (f[1] != 0.0f) ? ((f[0] != 0.0f) ? atan(f[1] / fabs(f[0])) : M_PI * f[1] / (2.0f * fabs(f[1]))) : 0.0f;
    if (f[0] < 0.0f) {
        a1 = ((f[1] < 0.0f) ? -1.0f : 1.0f) * M_PI - a1;
    }

    double a2 = (g[1] != 0.0f) ? ((g[0] != 0.0f) ? atan(g[1] / fabs(g[0])) : M_PI * g[1] / (2.0f * fabs(g[1]))) : 0.0f;
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

    printf("Computed shift: [%d, %d]\n", deltax, deltay);

    return 0;
}
