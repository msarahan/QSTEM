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

#include "pot_2d.hpp"

namespace QSTEM
{

C2DPotential::C2DPotential() : CPotential()
{
}

C2DPotential::C2DPotential(const ConfigReaderPtr &configReader) : CPotential(configReader)
{
	m_boxNz = 1;
}

void C2DPotential::Initialize()
{
}

void C2DPotential::Initialize(const ConfigReaderPtr &configReader)
{
}


void C2DPotential::DisplayParams()
{
  CPotential::DisplayParams();
  printf("* Potential calculation: 2D (real space method)");
}

void C2DPotential::AtomBoxLookUp(complex_tt &sum, int Znum, float_tt x, float_tt y, float_tt z, float_tt B) 
{
  float_tt dx, dy;
  int ix, iy;

  // does the atom box lookup or calculation
  CPotential::AtomBoxLookUp(sum, Znum, x, y, z, B);

  /***************************************************************
   * Do the trilinear interpolation
   */
  sum[0] = 0.0;
  sum[1] = 0.0;
  if (x*x+y*y+z*z > m_atomRadius2) {
    return;
  }
  x = fabs(x);
  y = fabs(y);
  ix = (int)(x/m_ddx);
  iy = (int)(y/m_ddy);
  dx = x-(float_tt)ix*m_ddx;
  dy = y-(float_tt)iy*m_ddy;
  
  if ((dx < 0) || (dy<0) ) {
    /* printf("Warning, dx(%g), dy(%g), dz(%g) < 0, (x=%g, y=%g, z=%g)\n",dx,dy,dz,x,y,z);
     */
    if (dx < 0) dx = 0.0;
    if (dy < 0) dy = 0.0;
  }
  
  if (m_atomBoxes[Znum]->B > 0) {
    sum[0] = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy][0]+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy][0])+
                       dy*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy+1][0]+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy+1][0]);
    sum[1] = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy][1]+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy][1])+
                       dy*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy+1][1]+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy+1][1]);
  }
  else {
    sum[0] = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->rpotential[0][ix][iy]+
                       dx*m_atomBoxes[Znum]->rpotential[0][ix+1][iy])+
                       dy*((1.0-dx)*m_atomBoxes[Znum]->rpotential[0][ix][iy+1]+
		       dx*m_atomBoxes[Znum]->rpotential[0][ix+1][iy+1]);
  }
}

bool C2DPotential::CheckAtomZInBounds(float_tt atomZ)
{
  /*
   * c = the thickness of the current slab.
   *
   * if the z-position of this atom is outside the potential slab
   * we won't consider it and skip to the next
   */
  return ((atomZ<m_c) && (atomZ>=0));
}


void C2DPotential::AddAtomToSlices(const atom &_atom, float_tt atomX, float_tt atomY,
                                               float_tt atomZ)
{
  // Note that if you override this method, you should do the following check to make sure the atom is in bounds.
  // skip atoms that are beyond the cell's boundaries
  if (!m_periodicZ)
    {
      if (atomZ > m_crystal->GetCZ() || atomZ<0) return;
    }

  // Calls parent class method, which in turn calls method below after computing ix, iy, iAtomZ
  AddAtomRealSpace(_atom, atomX, atomY, atomZ);
}

void C2DPotential::_AddAtomRealSpace(const atom &_atom, 
                                     float_tt atomBoxX, unsigned int ix, 
                                     float_tt atomBoxY, unsigned int iy, 
                                     float_tt atomZ, unsigned int iAtomZ)
{
  complex_tt dPot;
  // Coordinates in the total potential space
  unsigned iz;
  

  if (!m_periodicZ) {
    if (iAtomZ < 0) return;
    if (iAtomZ >= m_nslices) return;        
  }                
  iz = (iAtomZ+32*m_nslices) % m_nslices;         /* shift into the positive range */
  // x, y are the coordinates in the space of the atom box
  AtomBoxLookUp(dPot,_atom.Znum,atomBoxX,atomBoxY,0, m_tds ? 0 : _atom.dw);
  float_tt atomBoxZ = (double)(iAtomZ+1)*m_cz[0]-atomZ;

  unsigned idx=ix*m_ny+iy;

  /* split the atom if it is close to the top edge of the slice */
  if ((atomBoxZ<0.15*m_cz[0]) && (iz >0)) {
    m_trans[iz][idx] += 0.5*dPot;
    m_trans[iz-1][idx] += 0.5*dPot;
  }
  /* split the atom if it is close to the bottom edge of the slice */
  else {
    if ((atomBoxZ>0.85*m_cz[0]) && (iz < m_nslices-1)) {
      m_trans[iz][idx] += 0.5*dPot;
      m_trans[iz+1][idx] += 0.5*dPot;
    }
    else {
      m_trans[iz][idx] += dPot;
    }
  }
}

void C2DPotential::CenterAtomZ(const atom &_atom, float_tt &z)
{
  CPotential::CenterAtomZ(_atom, z);
  z += 0.5*m_sliceThickness;
}

} // end namespace QSTEM