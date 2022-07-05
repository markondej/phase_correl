# phase_correl

This is an example use of Phase Correlation algorithm to calculate shift between two images. Even though code has not been optimized for speed, this algorithm might be useful for real time image alignment. I successfully used it as a basis for my image processing researches in 2010 while working at Warsaw Univeristy of Technology.

Phase Correlation based image alignment was first described by C. D. Kuglin and D. C. Hines in the 1975 article "The Phase Correlation Image Alignment Method" (http://boutigny.free.fr/Astronomie/AstroSources/Kuglin-Hines.pdf).

## How it works

At first two monochromatic images are created, both containing squares of the same size but in different places. 

![Squares](squares.png)

Next, 2D Fourier transform is applied on each of them, normalised cross power spectrum is computed, and inverted 2D Fourier transform is applied. Finally algorightm searches for a peak in obtained matrix, where offset of peak found should match shift between images.

## Used FFT algorithm

For sake of simplicity used FFT is based on Radix-2 Cooley-Tukey alogrithm with input order calculation (not bit-reverse approach). Earlier versions (see: tags) have used recursive algorithm which (combined with large input images) might cause heap-overfill errors. Input order may be precalculated and more efficient radix may be used to improve performance.

## How to compile

Simply clone this repository on linux machine by typing:
```
git clone https://github.com/markondej/phase_correl.git
```
Then build executable binary with "make" command:
```
make
```
 
