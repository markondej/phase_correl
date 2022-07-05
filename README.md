# phase_correl

This is simple example of phase correlation use in order to compute shift between two images. Even though it has not been optimized for speed, this algorithm might be useful for real time image alignment. I used it as a base for my image processing researches in 2010 while working at Warsaw Univeristy of Technology.

Phase correlation feature matching was described by Yan H. and Liu J.G. in article "Roboust phase correlation based feature matching for image co-registration and DEM generation" (http://www.isprs.org/proceedings/XXXVII/congress/7_pdf/6_WG-VII-6/43.pdf).

## How it works

At first two monochromatic images are created, both containing squares of the same size but in different places. 

![Squares](squares.png)

Next, 2D Fourier transform is applied on each of them, normalised cross power spectrum is computed, and inverted 2D Fourier transform is applied. Finally algorightm searches for a peak in obtained matrix, where offset of peak found should match shift between images.

## Used FFT algorithm

For sake of simplicity used FFT is based on Radix-2 Cooley-Tukey alogrithm with input order calculation (not bit-reverse approach). Earlier versions (see: tags) have used recursive algorithm which (combined with large input images) might cause heap-overfill errors. Input order may be precalculated and more efficient radix may be used in order to improve performance.

## How to compile

Simply clone this repository on linux machine by typing:
```
git clone https://github.com/markondej/phase_correl.git
```
Then build executable binary with "make" command:
```
make
```
 
