# w_calculator
Calculator for angular correlation coefficients, for use in gamma ray spectroscopy. Written because there isn't an obvious resource online for quick calculation of Clebsch-Gordan, Racah, F- and A-coefficients (A2, A4 normalised to A0=1) relevant to gamma-gamma angular correlation measurements. Doesn't require anything other than the w_calculator.C file. Feel free to copy and edit your own version.


Copy-and-pasted from inside the .C file, explaining roughly what it contains:

// Calculates Clebsch-Gordan, Racah, F and A coefficients for gamma-gamma angular correlations. Also produces W-distribution between 0 and 90 degrees. Specify cascade by changing inside code (below)
//
// Compile and run with:
//  gcc -Wall -o w w_calculator.C -lm
//  ./w
//
// Produces .dat file containing plot of W(theta).
//
// t.berry@surrey.ac.uk
//
//  Multiple functions are written in order to calculate coefficients. The Wolfram definitions for the coefficients were used alongside the description of correlation coefficients in Taylor, Prato and Singh's 1971 paper:
//
//    Nuclear Data Tables A9, 1-83 (1971)
//
// This tabulation of A-coefficients is also useful as you can use it to check the validity of the result for given inputs. In its current form the code works perfectly so long as the inputs make physical sense. It won't necessarily prompt you if they are wrong, but it should be clear from the result whether the inputs have been entered incorrectly.
//
