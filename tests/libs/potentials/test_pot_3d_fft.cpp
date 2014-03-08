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
#define BOOST_TEST_MODULE Test3DFFTPotential
#include <boost/test/unit_test.hpp>
#include "boost_collection_close.hpp"

#include "potentials/pot_factory.hpp"

using namespace QSTEM;

struct PotFixture {
  PotFixture()
  {
    pot = CPotFactory::Get()->GetPotential(POTENTIAL3D, POTENTIALFFT);
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~PotFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  PotPtr pot;
};

BOOST_FIXTURE_TEST_SUITE(TestPotential3DFFT, PotFixture)

BOOST_AUTO_TEST_CASE(TestAtomBoxGeneration)
{
	// Calculation stores the potential as private member variable m_atomPot (map of atomic number to generated ComplexVector)
	//pot->GetAtomPotential3D();
}

BOOST_AUTO_TEST_SUITE_END()
