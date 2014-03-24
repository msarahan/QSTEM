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

#ifndef STRUCTURE_INTERFACE_H
#define STRUCTURE_INTERFACE_H

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include "stemtypes_fftw3.hpp"

#include <cstring>

namespace QSTEM
{

class IStructureReader;
typedef boost::shared_ptr<IStructureReader> StructureReaderPtr;
typedef StructureReaderPtr (*CreateStructureReaderFn)(const boost::filesystem::path &filename);

/**
Abstraction to read a crystal structure from a file.  Two things are essential:
- the matrix describing how the cell can be duplicated
- the list of atoms

This class is not meant to be used directly.  Instead, use the CCrystal class to manage your
crystal structures.  It will in turn use a reader appropriate to the filename you specify.
*/
class IStructureReader
{
public:
  virtual int ReadCellParams(RealVector &Mm)=0;
  virtual int ReadAtoms(std::vector<atom> &atoms)=0;
};


class IStructureWriter;
typedef boost::shared_ptr<IStructureWriter> StructureWriterPtr;
typedef StructureWriterPtr (*CreateStructureWriterFn)(const boost::filesystem::path &filename);

/**
Abstraction to write a crystal structure to a file.  Two things are essential:
- the matrix describing how the cell can be duplicated
- the list of atoms

This class is not meant to be used directly.  Instead, use the CCrystal class to manage your
crystal structures.  It will in turn use a writer appropriate to the filename you specify.
*/
class IStructureWriter
{
public:
  virtual int Write(std::vector<atom> &atoms, std::string run_id, float_tt ax, float_tt by, float_tt cz, 
	  float_tt alpha, float_tt beta, float_tt gamma)=0;
  virtual int Write(std::vector<atom> &atoms, const RealVector &Mm)=0;
};



}
#endif
