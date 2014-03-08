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
#define BOOST_TEST_MODULE TestStemExperiment
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

#include "experiments.hpp"

using namespace QSTEM;

struct ExperimentFixture {
  ExperimentFixture()
  {    
   configReader = CConfigReaderFactory::Get()->GetReader("stem_STO_4x4.qsc");
   expt = GetExperiment(configReader);
  }
  ~ExperimentFixture()
  { 
  }
  ConfigReaderPtr configReader;
  ExperimentPtr expt;
};

BOOST_FIXTURE_TEST_SUITE(TestExperimentSTEM, ExperimentFixture)

BOOST_AUTO_TEST_CASE( test_case )
{
  
}

BOOST_AUTO_TEST_SUITE_END()
