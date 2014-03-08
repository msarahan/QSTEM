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

#include "cfg.hpp"
#include "elTable.hpp"

#include "readparams.hpp"

namespace QSTEM
{

CCfgReader::CCfgReader(const boost::filesystem::path &filename)
{
  m_fp = fopen(filename.string().c_str(), "r" );
  if( m_fp == NULL ) {
    printf("Cannot open file %s\n",filename.string().c_str());
    m_isValid=false;
    return;
  }
  if (readparam(m_fp, ".NO_VELOCITY.",m_buf,1)) m_noVelocity = true; 
  else m_noVelocity = false;
  m_isValid=true;
}

CCfgReader::~CCfgReader()
{
  fclose(m_fp);
}

/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int CCfgReader::ReadCellParams(RealVector &Mm) {
  int i;
  double lengthScale;

  Mm.resize(9);

  resetParamFile(m_fp);  
  setComment('#');  
  if (readparam(m_fp, "A =",m_buf,1)) lengthScale=atof(m_buf.c_str());

  if (readparam(m_fp, "H0(1,1) =",m_buf,1)) Mm[0] = atof(m_buf.c_str());
  if (readparam(m_fp, "H0(1,2) =",m_buf,1)) Mm[1] = atof(m_buf.c_str());
  if (readparam(m_fp, "H0(1,3) =",m_buf,1)) Mm[2] = atof(m_buf.c_str());

  if (readparam(m_fp, "H0(2,1) =",m_buf,1)) Mm[3] = atof(m_buf.c_str());
  if (readparam(m_fp, "H0(2,2) =",m_buf,1)) Mm[4] = atof(m_buf.c_str());
  if (readparam(m_fp, "H0(2,3) =",m_buf,1)) Mm[5] = atof(m_buf.c_str());

  if (readparam(m_fp, "H0(3,1) =",m_buf,1)) Mm[6] = atof(m_buf.c_str());
  if (readparam(m_fp, "H0(3,2) =",m_buf,1)) Mm[7] = atof(m_buf.c_str());
  if (readparam(m_fp, "H0(3,3) =",m_buf,1)) Mm[8] = atof(m_buf.c_str());

  setComment('%');
  
  for (i=0;i<9;i++) Mm[i] *= lengthScale;

  return 0;
}

int CCfgReader::ReadAtoms(std::vector<atom> &atoms)
{
  unsigned ncoord;

  if (!m_isValid)
    throw std::runtime_error("Invalid structure file in CCfgReader.");

  resetParamFile(m_fp);
  if (readparam(m_fp, "Number of particles =",m_buf,1)) ncoord=atoi(m_buf.c_str());
  atoms.resize(ncoord);
  // This is important to be right here, because it sets the position of the file to be the first
  //     real atom entry.
  if (readparam(m_fp, "entry_count =",m_buf,1)) m_entryCount=atoi(m_buf.c_str());
  if (!m_noVelocity) m_entryCount+=3;
  m_atomData.resize(m_entryCount);
  for (unsigned i=0; i<ncoord; i++)
    {
      ReadNextAtom(&atoms[i]);
    }
  return 0;
}

/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .cfg file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* This function will return -1, if the end of file is reached prematurely, 0 otherwise.
*******************************************************************************/

int CCfgReader::ReadNextAtom(atom *newAtom) {
  char *str = NULL;
  char buf[NCMAX];

  if (m_fp == NULL) {
    printf("Invalid CFG file!\n");
    return -1;
  }

  if (fgets(buf,NCMAX,m_fp) == NULL) return -1;
  /* check, if this is a new mass number */
  str = strnext(buf," \t");
  if ((atof(buf) >= 1.0) && ((str==NULL) || (*str == '#'))) {
    m_lastReadMass = atof(buf);
    // printf("nV: %d, eC: %d (%g)\n",noVelocityFlag, entryCount,atof(m_buf));
    if (fgets(buf,NCMAX,m_fp) == NULL) return -1;    
    m_lastReadElement = getZNumber(buf); 
    // printf("*** found element %d (%s %d) ***\n",element,m_buf,strlen(m_buf));
    if (fgets(buf,NCMAX,m_fp) == NULL) return -1;
  }
  str = buf;
  // skip leading spaces:
  while (strchr(" \t",*str) != NULL) str++; 
  for (unsigned j=0;j<m_entryCount;j++) {
    if (str==NULL) {
      printf("readNextCFGatom: Error: incomplete data line: >%s<\n",buf);
      return -1;
    }
    m_atomData[j] = atof(str); str=strnext(str," \t");
  }

  newAtom->x    = m_atomData[0];
  newAtom->y    = m_atomData[1];
  newAtom->z    = m_atomData[2];
  // Temporary initialization: dw, occ, q are read if they exist
  newAtom->dw   = 0.45*28.0/m_lastReadMass;	
  newAtom->occ  = 1.0;
  newAtom->q    = 0.0;

  newAtom->Znum = m_lastReadElement;
  newAtom->mass = m_lastReadMass;
  // read the DW-factor
  if (m_entryCount > 3+3*(1-(int)m_noVelocity)) 
    newAtom->dw = m_atomData[3+3*(1-(int)m_noVelocity)];
  // read the atom's occupancy:
  if (m_entryCount > 4+3*(1-(int)m_noVelocity)) 
    newAtom->occ = m_atomData[4+3*(1-(int)m_noVelocity)];
  // read the atom's charge:
  if (m_entryCount > 5+3*(1-(int)m_noVelocity)) 
    newAtom->q = m_atomData[5+3*(1-(int)m_noVelocity)];
  // printf("Atom: %d (%g %g %g), occ=%g, q=%g\n",newAtom->Znum,newAtom->x,newAtom->y,newAtom->z,newAtom->occ,newAtom->q);	

  return 0;
}

CCfgWriter::CCfgWriter(const boost::filesystem::path &path)
  : m_basePath(path)
{
}

CCfgWriter::~CCfgWriter()
{
}

int CCfgWriter::Write(std::vector<atom> &atoms, std::string run_id, float_tt ax, float_tt by, float_tt cz, 
	  float_tt alpha, float_tt beta, float_tt gamma) {
  int j;
  std::string elem;

  if (atoms.size() < 1) {
    printf("Atom array empty - no file written\n");
    return 1;
  }
  
  std::string outputFile = m_basePath.stem().string();
  if (run_id!="")
    {
      outputFile+="_";
      outputFile+=run_id;
    }
  outputFile+=m_basePath.extension().string();
  FILE *fp = fopen(outputFile.c_str(), "w" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",outputFile.c_str());
  }

  fprintf(fp,"Number of particles = %d\n",atoms.size());
  fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
  fprintf(fp,"H0(1,1) = %g A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n",ax);
  fprintf(fp,"H0(2,1) = 0 A\nH0(2,2) = %g A\nH0(2,3) = 0 A\n",by);
  fprintf(fp,"H0(3,1) = 0 A\nH0(3,2) = 0 A\nH0(3,3) = %g A\n",cz);
  fprintf(fp,".NO_VELOCITY.\nentry_count = 6\n");

  //printf("ax: %g, by: %g, cz: %g n: %d\n",m_ax,m_by,m_cz,atoms.size());


  elem = getSymbol(atoms[0].Znum);
  fprintf(fp,"%g\n%s\n",atoms[0].mass,elem.c_str());
  fprintf(fp,"%g %g %g %g %.4f %.4f\n",atoms[0].x,atoms[0].y,atoms[0].z,
          atoms[0].dw,atoms[0].occ,atoms[0].q);

  for (j=1;j<atoms.size();j++) {
    if (atoms[j].Znum != atoms[j-1].Znum) {
		elem = getSymbol(atoms[j].Znum);
      fprintf(fp,"%g\n%s\n",atoms[j].mass,elem.c_str());
      // printf("%d: %g\n%s\n",j,2.0*atoms[j].Znum,elem);
    }
    fprintf(fp,"%g %g %g %g %.4f %.4f\n",atoms[j].x,atoms[j].y,atoms[j].z,
            atoms[j].dw,atoms[j].occ,atoms[j].q);
    // if (atoms[j].occ != 1) printf("Atom %d: occ = %g\n",j,atoms[j].occ);
  }
  fclose(fp);

  return 1;
}


#define MIN_EDGE_LENGTH 5.18 /* minimal allowed edge length in A
* going below this limit will crash 
* AtomEye.
*/

} // end namespace QSTEM
