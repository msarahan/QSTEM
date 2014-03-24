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

#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

// Inlcude this so that we know what an ExperimentPtr is
#include "experiments/experiment_interface.hpp"
// Include this because we use ConfigReaderPtr below.
#include "config_IO/config_reader_factory.hpp"

namespace QSTEM
{

ExperimentPtr QSTEM_HELPER_DLL_EXPORT GetExperiment(ConfigReaderPtr &configReader);

}

#endif
