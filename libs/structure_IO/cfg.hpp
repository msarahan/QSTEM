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

#ifndef CFG_STRUCTURE_H
#define CFG_STRUCTURE_H

#include "structureInterface.hpp"
#include <boost/filesystem.hpp>
namespace QSTEM
{
class QSTEM_HELPER_DLL_EXPORT CCfgReader : public IStructureReader
{
public:
  CCfgReader(const boost::filesystem::path &structure_file);
  ~CCfgReader();
  
  int ReadCellParams(RealVector &Mm);
  int ReadAtoms(std::vector<atom> &atoms);
protected:
  int ReadNextAtom(atom *newAtomcom);
private:
  FILE *m_fp;
  std::vector<float_tt> m_atomData; // One line of atom data read from a file
  bool m_noVelocity;     // If true, the file does not provide velocity information
  unsigned m_entryCount; // The number of columns on any given line of data from the file
  bool m_isValid;  // If true, this object has a valid file to read atoms from.
  std::string m_buf;
  float_tt m_lastReadMass;
  unsigned m_lastReadElement;
private:
  friend class CStructureReaderFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static StructureReaderPtr Create(const boost::filesystem::path &file){
    return StructureReaderPtr(new CCfgReader(file));}
};

class CCfgWriter : public IStructureWriter
{
public:
  CCfgWriter(const boost::filesystem::path &file);
  ~CCfgWriter();

  int Write(std::vector <atom> &atoms, std::string run_id, float_tt ax, float_tt by, float_tt cz, 
	  float_tt alpha, float_tt beta, float_tt gamma);
  int Write(std::vector<atom> &atoms, const RealVector &Mm);
private:
  boost::filesystem::path m_basePath;
  FILE *OpenFile(const std::string &run_id);
  void WriteAtoms(FILE *fp, std::vector<atom> atoms);
private:
  friend class CStructureWriterFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static StructureWriterPtr Create(const boost::filesystem::path &file){
    return StructureWriterPtr(new CCfgWriter(file));
  }  
};

}
#endif

















