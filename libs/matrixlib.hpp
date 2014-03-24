/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
	Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MATRIXLIB_H
#define MATRIXLIB_H

#include "stemtypes_fftw3.hpp"

#define PI 3.1415926535898
#define PI180 1.7453292519943e-2

namespace QSTEM
{

void ludcmp(RealVector &a, int n, int *indx, RealVector &d);
void lubksb(RealVector &a, int n, int *indx, float_tt b[]);
float_tt det_3x3 (const RealVector &mat);
void inverse_3x3 (RealVector &res, const RealVector &a);
void trans_3x3 (RealVector &Mt, const RealVector &Ms);

// svdcmp1 uses the NR unit-offset vectors :-(
void svdcmp1(RealVector &a, int m, int n, float_tt w[], RealVector &v);
float_tt pythag(float_tt a, float_tt b);

/* vector functions:
 */
void crossProduct(const RealVector &a, const RealVector &b, RealVector &c);
float_tt dotProduct(const RealVector &a, const RealVector &b);
float_tt findLambda(plane *p, RealVector &point, int revFlag);
void showMatrix(RealVector &M,int Nx, int Ny,char *name);
void vectDiff(RealVector &a, float_tt *b, RealVector &c,int revFlag);
float_tt vectLength(RealVector &vect);
void makeCellVect(grainBox *grain, RealVector &vax, RealVector &vby, RealVector &vcz);
//void makeCellVectMuls(MULS *muls, RealVector &vax, RealVector &vby, RealVector &vcz);
void rotateVect(const RealVector &vectIn,RealVector &vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
void rotateMatrix(const RealVector &matrixIn,RealVector &matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
void populateRotationMatrix(RealVector &Mm, float_tt phi_x, float_tt phi_y, float_tt phi_z);

/* |vect| */
float_tt vectLength(RealVector &vect);

template <typename T>
void rotateVect(const T &vectIn, T &vectOut, const RealVector &Mm)
{
	matrixProduct(Mm, 3, 3, vectIn, 3, 1, vectOut);
}

  /* c = a*b */
  template <typename T>
  void matrixProduct(const T &a,int Nxa, int Nya, const T &b,int Nxb, int Nyb, T &c) {
    int i,j,k;

    if (Nya != Nxb) {
      printf("multiplyMatrix: Inner Matrix dimensions do not agree!\n");
      return;
    }

    for (i=0;i<Nxa;i++) for (j=0;j<Nyb;j++) {
      c[i*Nyb+j] = 0.0;
      for (k=0;k<Nya;k++) c[i*Nyb+j] += a[i*Nya+k] * b[k*Nyb+j];
    }
}


} // end namespace QSTEM

#endif
