#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IMAGE_WIDTH 256
#define IMAGE_HEIGHT 128

void ditfft2(double *input, double *output, int size, int s, int invert) {
    double temp[4], fi, kfi, cosfi, sinfi;
    int k, k2;

    if (size > 1) {
        ditfft2(&input[0], &output[0], size >> 1, s << 1, invert);
        ditfft2(&input[s << 1], &output[size], size >> 1, s << 1, invert);

        if (!invert) {
            fi = -2.0 * M_PI / (double)size; kfi = 0.0;
        } else {
            fi = 2.0 * M_PI / (double)size; kfi = 0.0;
        }
        for (k = 0; k < size >> 1; k++) {
            cosfi = cos(kfi); sinfi = sin(kfi); k2 = k << 1;
            temp[0] = output[k2]; temp[1] = output[k2 + 1];
            temp[2] = output[size + k2]; temp[3] = output[size + k2 + 1];

            output[k2] = temp[0] + temp[2] * cosfi - temp[3] * sinfi;
            output[k2 + 1] = temp[1] + temp[2] * sinfi + temp[3] * cosfi;
            output[size + k2] = temp[0] - temp[2] * cosfi + temp[3] * sinfi;
            output[size + k2 + 1] = temp[1] - temp[2] * sinfi - temp[3] * cosfi;
            kfi += fi;
        }
    }
    else {
        output[0] = input[0];
        output[1] = input[1];
    }
}

void fft(double *input, double *output, int size) {
    ditfft2(input, output, size, 1, 0);
}

void ifft(double *input, double *output, int size) {
    ditfft2(input, output, size, 1, 1);
}

void fft2D(double *input, int width, int height) {
    int i, j;

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
    int i, j;

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

void computeShift(unsigned char *image1, unsigned char *image2, int width, int height, int *tx, int *ty) {
    int i, j;

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
        double a = fft_input1[i << 1];
        double b = fft_input1[(i << 1) + 1];
        double c = fft_input2[i << 1];
        double d = fft_input2[(i << 1) + 1];

        double a1, a2;

        if (b != 0.0) {
            if (a != 0.0) {
                a1 = atan(fabs(b / a));
            }
            else {
                a1 = M_PI / 2.0;
            }
        }
        else a1 = 0.0;

        if ((a < 0.0) && (b > 0.0)) a1 = M_PI - a1;
        if ((a < 0.0) && (b < 0.0)) a1 = M_PI + a1;
        if ((a > 0.0) && (b < 0.0)) a1 = 2.0 * M_PI - a1;

        if (d != 0.0) {
            if (c != 0.0) {
                a2 = atan(fabs(d / c));
            }
            else {
                a2 = M_PI / 2.0;
            }
        }
        else a2 = 0.0;

        if ((c < 0.0) && (d > 0.0)) a2 = M_PI - a2;
        if ((c < 0.0) && (d < 0.0)) a2 = M_PI + a2;
        if ((c > 0.0) && (d < 0.0)) a2 = 2.0 * M_PI - a2;

        fft_output[i << 1] = cos(a1 - a2);
        fft_output[(i << 1) + 1] = sin(a1 - a2);
    }

    // Perform inversed 2D FFT on obtained matrix
    ifft2D(fft_output, width, height);

    // Search for peak
    double max = 0.0; *tx = 0; *ty = 0;
    for (j = 0; j < height; j++)
        for (i = 0; i < width; i++) {
            double d = sqrt(pow(fft_output[(i + j * width) << 1], 2) + pow(fft_output[((i + j * width) << 1) + 1], 2));
            if (d > max) {
                max = d; *tx = i; *ty = j;
            }
        }

    if (*tx > width >> 1)
        *tx = *tx - width;
    if (*ty > height >> 1)
        *ty = *ty - height;

    free(fft_input1);
    free(fft_input2);
    free(fft_output);
}

int main()
{
    unsigned char image1[IMAGE_WIDTH * IMAGE_HEIGHT];
    unsigned char image2[IMAGE_WIDTH * IMAGE_HEIGHT];

    int i, j, tx, ty;

    // Generate pair of images
    for (j = 0; j < IMAGE_HEIGHT; j++)
        for (i = 0; i < IMAGE_WIDTH; i++) {
            if ((i > 16) && (i < 80) && (j > 32) && (j < 96))
                image1[i + j * IMAGE_WIDTH] = 128;
            else
                image1[i + j * IMAGE_WIDTH] = 0;

            if ((i > 8) && (i < 72) && (j > 40) && (j < 104))
                image2[i + j * IMAGE_WIDTH] = 16;
            else
                image2[i + j * IMAGE_WIDTH] = 255;
        }

    computeShift(image1, image2, IMAGE_WIDTH, IMAGE_HEIGHT, &tx, &ty);

    printf("Computed shift: [%d, %d]\n", tx, ty);

    return 0;
}
