#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

void radix2fft(const double *input, double *output, size_t stride) {
    output[0] = input[0] + input[stride << 1];
    output[1] = input[1] + input[(stride << 1) + 1];

    output[2] = input[0] - input[(stride << 1)];
    output[3] = input[1] - input[(stride << 1) + 1];
}

void sum2fft(const double *input, double *output, size_t size, int inverse) {
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
        output[(k + size) << 1] = temp[0] - cosfi * temp[2] + sinfi * temp[3];
        output[((k + size) << 1) + 1] = temp[1] - sinfi * temp[2] - cosfi * temp[3];
        kfi += dfi;
    }
}

size_t *gen2fftorder(size_t size) {
    size_t *offset = (size_t *)malloc(sizeof(size_t) * size),
        stride = 1, step;

    memset(offset, 0xff, sizeof(size_t) * size);
    offset[0] = 0;

    for (step = size >> 1; step > 0; step >>= 1) {
        unsigned base = 0, i;
        for (i = 0; i < size; i += step)
            if (offset[i] != SIZE_MAX)
                base = offset[i];
            else
                offset[i] = base + stride;
        stride <<= 1;
    }

    return offset;
}

void dit2fft(const double *input, double *output, size_t size, int inverse) {
    size_t *offset = gen2fftorder(size >> 1),
        stride = 1, i;

    while (size > 2) {
        stride <<= 1;
        size >>= 1;
    }

    for (i = 0; i < stride; i++)
        radix2fft(&input[offset[i] << 1], &output[(i * 2) << 1], stride);

    free(offset);

    while (stride > 1) {
        stride >>= 1;
        size <<= 1;
        for (i = 0; i < stride; i++)
            sum2fft(&output[(i * size) << 1], &output[(i * size) << 1], size >> 1, inverse);
    }

    if (inverse)
        for (i = 0; i < size; i++) {
            output[i << 1] /= (double)size;
            output[(i << 1) + 1] /= (double)size;
        }
}

int fft2D(double *data, unsigned long width, unsigned long height, int inverse) {
    unsigned long i, j;

    double *fft_input = (double *)malloc((sizeof(double) * width) << 1),
        *fft_output = (double *)malloc((sizeof(double) * width) << 1);

    if (!fft_input || !fft_output)
        return -1;

    if (inverse)
        goto horizontal_fft;

vertical_fft:
    for (j = 0; j < height; j++) {
        size_t offset = j * width;
        for (i = 0; i < width; i++) {
            fft_input[i << 1] = data[offset << 1];
            fft_input[(i << 1) + 1] = data[(offset << 1) + 1];
            offset++;
        }
        dit2fft(fft_input, fft_output, width, inverse);
        offset = j * width;
        for (i = 0; i < width; i++) {
            data[offset << 1] = fft_output[i << 1];
            data[(offset << 1) + 1] = fft_output[(i << 1) + 1];
            offset++;
        }
    }

    if (inverse)
        goto done;

horizontal_fft:
    for (i = 0; i < width; i++) {
        size_t offset = i;
        for (j = 0; j < height; j++) {
            fft_input[j << 1] = data[offset << 1];
            fft_input[(j << 1) + 1] = data[(offset << 1) + 1];
            offset += width;
        }
        dit2fft(fft_input, fft_output, height, inverse);
        offset = i;
        for (j = 0; j < height; j++) {
            data[offset << 1] = fft_output[j << 1];
            data[(offset << 1) + 1] = fft_output[(j << 1) + 1];
            offset += width;
        }
    }

    if (inverse)
        goto vertical_fft;

done:
    free(fft_input);
    free(fft_output);

    return 0;
}

void computeNormalized(const double *f, const double *g, double *r) {
    double a1 = (f[1] != 0.0) ? ((f[0] != 0.0) ? atan(f[1] / f[0]) : copysign(M_PI, f[1]) / 2.0) : 0.0;
    if (f[0] < 0.0)
        a1 = copysign(M_PI, f[1]) + a1;

    double a2 = (g[1] != 0.0) ? ((g[0] != 0.0) ? atan(g[1] / g[0]) : copysign(M_PI, g[1]) / 2.0) : 0.0;
    if (g[0] < 0.0)
        a2 = copysign(M_PI, g[1]) + a2;

    r[0] = cos(a1 - a2);
    r[1] = sin(a1 - a2);
}

int computeShift(const unsigned char *image1, const unsigned char *image2, unsigned long width, unsigned long height, long *deltax, long *deltay) {
    if (!width || (width & (width - 1)) || !height || (height & (height - 1)))
        return -1;

    unsigned long i, j;
    int ret = 0;

    double *fft_input1 = (double *)malloc((sizeof(double) * width * height) << 1),
        *fft_input2 = (double *)malloc((sizeof(double) * width * height) << 1),
        *fft_output = (double *)malloc((sizeof(double) * width * height) << 1);

    if (!fft_input1 || !fft_input2 || !fft_output)
        return -1;

    // Convert image pixels to complex number format, use only real part
    for (i = 0; i < width * height; i++) {
        fft_input1[i << 1] = (double)image1[i];
        fft_input2[i << 1] = (double)image2[i];

        fft_input1[(i << 1) + 1] = 0.0;
        fft_input2[(i << 1) + 1] = 0.0;
    }

    // Perform 2D FFT on each image
    ret = fft2D(fft_input1, width, height, 0);
    if (ret)
        goto clean;
    ret = fft2D(fft_input2, width, height, 0);
    if (ret)
        goto clean;

    // Compute normalized cross power spectrum
    for (i = 0; i < width * height; i++)
        computeNormalized(&fft_input1[i << 1], &fft_input2[i << 1], &fft_output[i << 1]);

    // Perform inverse 2D FFT on obtained matrix
    ret = fft2D(fft_output, width, height, 1);
    if (ret)
        goto clean;

    // Search for peak
    size_t offset = 0;
    double max = 0.0; *deltax = 0; *deltay = 0;
    for (j = 0; j < height; j++)
        for (i = 0; i < width; i++) {
            double d = sqrt(pow(fft_output[offset << 1], 2) + pow(fft_output[(offset << 1) + 1], 2));
            if (d > max) {
                max = d;
                *deltax = i;
                *deltay = j;
            }
            offset++;
        }

    if (*deltax > width >> 1)
        *deltax -= width;
    if (*deltay > height >> 1)
        *deltay -= height;

clean:
    free(fft_input1);
    free(fft_input2);
    free(fft_output);

    return ret;
}

int main()
{
    unsigned char image1[256 * 128];
    unsigned char image2[256 * 128];

    unsigned long i, j;
    long deltax, deltay;

    // Generate pair of images
    for (j = 0; j < 128; j++)
        for (i = 0; i < 256; i++) {
            unsigned offset = i + j * 256;
            if ((i >= 16) && (i < 76) && (j >= 32) && (j < 92))
                image1[offset] = 128;
            else
                image1[offset] = 0;

            if ((i >= 8) && (i < 68) && (j >= 40) && (j < 100))
                image2[offset] = 16;
            else
                image2[offset] = 255;
        }

    int ret = computeShift(image1, image2, 256, 128, &deltax, &deltay);
    if (ret) {
        fprintf(stderr, "Operation failed");
        return EXIT_FAILURE;
    }

    printf("Calculated shift: [%ld, %ld]\n", deltax, deltay);
    return EXIT_SUCCESS;
}
