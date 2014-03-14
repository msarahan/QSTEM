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

#include "stemtypes_fftw3.hpp"

#ifndef EXPERIMENT_INTERFACE_H
#define EXPERIMENT_INTERFACE_H

namespace QSTEM
{

/**
Experiments are meant to be outlines and assemblies of simulation components.  They
tie together lesser parts, and take care of how the beam starts and where/how it gets detected.
For example, a multislice experiment will start by creating Wavefunction and Crystal objects.  
From there, it will compute the Potential using the Crystal.  It will take care of propagating the Wavefunction
through the Potential, and finally save output.

Experiments also encapsulate state such as probe position and any averaging.
*/
class IExperiment
{
public:
  virtual void DisplayParams(){}; 
  virtual void CheckParams()=0;
  virtual void DescribeParams()=0;
  virtual void Run()=0;
protected:
  unsigned m_printLevel;
  unsigned m_saveLevel;
};

typedef boost::shared_ptr<IExperiment> ExperimentPtr;

}

#endif
