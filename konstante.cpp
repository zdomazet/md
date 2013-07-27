#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <valarray>
#include <fstream>
#include <complex>

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

#define u 1.660538e-3   //ng
#define pi 3.1415926
#define kb 1.3806503e-5  //nm^2 ng ps^-2 K^-1
const double PI=M_PI;

double m=u*39.948;    //ng
double T=100.;        //K
double dt=0.00001;       //ps
double a0=0.5;         //nm

int n=100000;           //broj vremenskih koraka 

double sigma=0.3345;  //nm
double epsilon=125.7*kb;

double A=4*epsilon*pow(sigma,12);
double B=4*epsilon*pow(sigma,6);

double cutoff=4.*sigma;

typedef struct {
        double r[3];
        double v[3];
        double F[3];
        } cestica;

int size=6;
int N=4*(int)pow(size,3);
double V=pow((size+1)*a0, 3);

double rdf_dr=0.001;
int n_rdf=(int)2*(size+1)*a0/rdf_dr;

double U=1.5*N*kb*T;               //termalna energija
double U_trunc=A/pow(cutoff, 12)-B/pow(cutoff,6);
