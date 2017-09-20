
//------------------------------------------------------------------------
// Calculates Clebsch-Gordan, Racah, F and A coefficients
// for gamma-gamma angular correlations.
// Also produces W-distribution between 0 and 90 degrees.
//
// Specify cascade by changing inside code (below)
//
// Compile and run with:
//  gcc -Wall -o w w_calculator.C -lm
//  ./w
//
// Produces .dat file containing plot of W(theta).
//
// t.berry@surrey.ac.uk
// -----------------------------------------------------------------------
//
//  Multiple functions are written in order to calculate coefficients.
//  The Wolfram definitions for the coefficients were used alongside the
//  description of correlation coefficients in Taylor, Prato and Singh's
//  1971 paper:
//
//    Nuclear Data Tables A9, 1-83 (1971)
//
//  This tabulation of A-coefficients is also useful as you can use it to
//  check the validity of the result for given inputs. In its current form
//  the code works perfectly so long as the inputs make physical sense.
//  It won't necessarily prompt you if they are wrong, but it should be
//  clear from the result whether the inputs have been entered incorrectly.
//
//------------------------------------------------------------------------



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float fac(float input);
float tsum6(float j1, float j2, float j3, float J1, float J2, float J3);
float tsum3(float j1, float j2, float j, float m1, float m2, float m);
double delta(float a, float b, float c);
double clebsch(float j1, float j2, float j,  float m1, float m2, float m);
double racah(float j1, float j2, float j3, float J1, float J2, float J3);
double fcalc(float k, float La, float Lb, float ja, float jb);
float max(float a1, float a2, float a3, float a4);
float min(float b1, float b2, float b3);

int main()
{

  float I0=0.0, I1=0.0, I2=0.0;
  int MP1=0, MP2=0;
  float L1=0.0, L1p=0.0, L2=0.0, L2p=0.0;
  float d1=0.0, d2=0.0;

// Cascade is 'visualised' below.
// Change level spins (I0,I1,I2), transition multipolarities (MP1,MP2) and mixing ratios (d1,d2).
// Currently only L+1 mixing (L1p,L2p) is accounted for, but feel free to change that manually below.
//
//            CASCADE
//
              I0=4.0;
//          -----------
//               |
//               |
//               |
//               |
              MP1=2;
              d1=0.0;
//               |
//               |
//               |
//              \ /
//               v
              I1=2.0;
//          -----------
//               |
//               |
//               |
//               |
//              \ /
//               v
              MP2=2;
              d2=0.0;
//               |
//               |
//               |
//              \ /
//               v
              I2=0.0;
//          -----------


  L1 = (float) MP1;
  L1p = L1+1.0;
  L2 = (float) MP2;
  L2p = L2+1.0;

  // Printing inputs for clarity.
  printf("\n\n Calculation of A2,A4 coefficients for angular correlations\n\n");
  printf("I0 %f  L1 %f  L1p %f  d1 %f  I1 %f  L2 %f  L2p %f  d2 %f  I2 %f \n",I0,L1,L1p,d1,I1,L2,L2p,d2,I2);

  double F2[6];
  double F4[6];

  double A2;
  double A4;



  // Calculation of F-coefficients from inputs above.
  F2[0] = fcalc(2.0,L1 ,L1 ,I0,I1);
  F2[1] = fcalc(2.0,L1 ,L1p,I0,I1);
  F2[2] = fcalc(2.0,L1p,L1p,I0,I1);
  printf("\n");
  F4[0] = fcalc(4.0,L1 ,L1 ,I0,I1);
  F4[1] = fcalc(4.0,L1 ,L1p,I0,I1);
  F4[2] = fcalc(4.0,L1p,L1p,I0,I1);
  printf("\n");
  F2[3] = fcalc(2.0,L2 ,L2 ,I2,I1);
  F2[4] = fcalc(2.0,L2 ,L2p,I2,I1);
  F2[5] = fcalc(2.0,L2p,L2p,I2,I1);
  printf("\n");
  F4[3] = fcalc(4.0,L2 ,L2 ,I2,I1);
  F4[4] = fcalc(4.0,L2 ,L2p,I2,I1);
  F4[5] = fcalc(4.0,L2p,L2p,I2,I1);
  printf("\n");

  // Main calculation of theoretical A2,A4 coefficients.

  A2 = ((1/(1+d1*d1))*(F2[0]+(pow(-1,L1-L1p)*2*d1*F2[1])+(d1*d1*F2[2]))) * ((1/(1+d2*d2))*(F2[3]+(2*d2*F2[4])+(d2*d2*F2[5])));

  A4 = ((1/(1+d1*d1))*(F4[0]+(pow(-1,L1-L1p)*2*d1*F4[1])+(d1*d1*F4[2]))) * ((1/(1+d2*d2))*(F4[3]+(2*d2*F4[4])+(d2*d2*F4[5])));

  printf("\n F2 values: %lf %lf %lf %lf %lf %lf \n",F2[0],F2[1],F2[2],F2[3],F2[4],F2[5]);
  printf("\n F4 values: %lf %lf %lf %lf %lf %lf \n",F4[0],F4[1],F4[2],F4[3],F4[4],F4[5]);

  printf("\n A2 result is %lf     A4 result is %lf\n\n",A2,A4);

  // Prints the value of W between theta=0-90 degrees
  char *output_filename = (char *)malloc(40*sizeof(char));
  sprintf(output_filename,"wcalc_%.1f_%.1f_%.1f.dat",I0,I1,I2);
  FILE * outfile = fopen(output_filename,"w");

  double W=0.0, P2=0.0, P4=0.0, x=0.0;

  for(int theta = 0; theta <= 90; theta++)
  {
    x = cos(theta*M_PI/180)*cos(theta*M_PI/180);

    P2 = 0.5000 * (3.0*x - 1.0);
    P4 = 0.1250 * (35.0*x*x - 30.0*x + 3.0);

    W = 1.0 + A2*P2 + A4*P4;

    fprintf(outfile, "%d   %lf\n",theta, W);
  }

  fclose(outfile);
  return(0);

}

double fcalc(float k, float La, float Lb, float ja, float jb)
{
  double f;

  // F is calculated from Clebsch-Gordan and Racah coefficients, calculated themselves in a separate external function.
  // The values of F are printed here mainly to check it's working properly.
  printf("\n F%d(%f %f %f %f) \n",(int)k,La,Lb,ja,jb);
  f = pow(-1.0,(ja-jb-1.0)) * sqrt((1.0+2.0*La)*(1.0+2.0*Lb)*(1.0+2.0*jb)) * clebsch(La,Lb,k,1.0,-1.0,0.0) * racah(jb,jb,k,La,Lb,ja);

  return(f);
}


double clebsch(float j1, float j2, float j, float m1, float m2, float m)
{
  double coeff;
  float s = j+j1+j2;

  coeff = sqrt( (1.0+2.0*j)*fac(s-2.0*j)*fac(s-2.0*j1)*fac(s-2.0*j2)*fac(j1+m1)*fac(j1-m1)*fac(j2+m2)*fac(j2-m2)*fac(j+m)*fac(j-m)/fac(s+1.0) );

  coeff *= tsum3(j1,j2,j,m1,m2,m);

  printf(" clebsch %lf ",coeff);

  return(coeff);
}


double racah(float j1, float j2, float j3, float J1, float J2, float J3)
{
  double coeff;

  coeff = sqrt(delta(j1, j2, j3)*delta(j1,J1,J3)*delta(j2,J2,J3)*delta(J1,J2,j3)) * tsum6(j1,j2,j3,J1,J2,J3);

  printf(" racah %lf ",coeff);

  return(coeff);
}

// Delta function is require for Racah coefficient calculation.
double delta(float a, float b, float c)
{
  double delta;

  delta = fac(a+b-c) * fac(a-b+c) * fac(b+c-a) / fac(a+b+c+1);

  return(delta);
}

// Factorials are required for Clebsch-Gordan coefficients.
float fac(float input)
{
  float factorial=1.0;

  if(input==0.0) factorial=1.0;

  else
  {
    for(int i = (int) input; i>0; i-=1)
    {
      factorial *= i;
    }
  }

  return(factorial);

}

// Sum required for Racah calculation.
float tsum6(float j1, float j2, float j3, float J1, float J2, float J3)
{

  float z = 0.0;
  float tlower=0.0, tupper=0.0;
  float tsum = 0.0;
  float f = 0.0;


  float a1=j1+j2+j3;
  float a2=J1+J2+j3;
  float a3=j1+J1+J3;
  float a4=j2+J2+J3;

  float b1=j1+j2+J1+J2;
  float b2=j1+J2+j3+J3;
  float b3=j2+J1+j3+J3;

  // Definition of max, min t-values to shorten calculation. See Wolfram site for full explanation.
  tlower = max(a1,a2,a3,a4);

  tupper = min(b1,b2,b3);

  for(z=tlower; z<=tupper; z+=1.0)
  {
    if((a1<=z)&&(a2<=z)&&(a3<=z)&&(a4<=z)&&(b1>=z)&&(b2>=z)&&(b3>=z))
    {
      f = fac(z+1.0) / ( fac(z-a1) * fac(z-a2) * fac(z-a3) * fac(z-a4) * fac(b1-z) * fac(b2-z) * fac(b3-z) );

      tsum += (pow(-1,(z+b1)) * f);

    }

  }

  return(tsum);
}

// Sum required for Clebsch-Gordan calculation.
float tsum3(float j1, float j2, float j, float m1, float m2, float m)
{

  float t = 0.0;
  float tsum = 0.0;
  float f = 0.0;
  int tint=0;

  // The value of 50 is arbitrary here; just needs to be high enough to cover all possible combinations easily.
  for(t=0.0; t<50.0; t+=1.0)
  {
    if(((j-j2+t+m1)>=0.0)&&((j-j1+t-m2)>=0.0)&&((j1+j2-j-t)>=0.0)&&((j1-t-m1)>=0.0)&&((j2-t+m2)>=0.0))
    {
      f = fac(t) * fac(j-j2+t+m1) * fac(j-j1+t-m2) * fac(j1+j2-j-t) * fac(j1-t-m1) * fac(j2-t+m2);

      tint = (int) t;
      tsum += pow(-1,tint) / f;
    }
  }

  return(tsum);
}

// Maximum and minimum is required in tsum6 to save time on loop. Defined in their own functions below.
float max(float a1, float a2, float a3, float a4)
{
  float result=0.0;

  result = 0.5*((a1+a2)+abs(a1-a2));
  result = 0.5*((result+a3)+abs(result-a3));
  result = 0.5*((result+a4)+abs(result-a4));

  return(result);
}

float min(float b1, float b2, float b3)
{
  float result=0.0;

  result = 0.5*((b1+b2)-abs(b1-b2));
  result = 0.5*((result+b3)-abs(result-b3));

  return(result);
}

