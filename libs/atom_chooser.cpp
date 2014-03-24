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

#include "atom_chooser.hpp"

using namespace QSTEM;

AtomChooser::AtomChooser(std::vector<atom> siteAtomList)
	: m_atoms(siteAtomList)
{
	RealVector probabilities = GetOccupancyVector();
	// sum occupancies
	float_tt occupancy_sum=SumOccupancy();
	// There are 3 possible outcomes:
	// occupancies sum to 1: no change (no vacancies)
	if (fabs(occupancy_sum-1)<1E-3)
	{
		// do nothing - the vector is fine as is
	}
	// occupancies sum to greater than 1: probabilities normalized to range from 0 to 1 (no vacancies)
	else if (occupancy_sum>1)
	{
		// normalize such that sum is 1 (divide each by the sum)
		RealVector::iterator prob=probabilities.begin(), end=probabilities.end();
		for (prob; prob!=end; ++prob)
		{
			(*prob)/=occupancy_sum;
		}
	}
	// occupancies sum to less than one: remainder is probability of vacancy
	else
	{
		atom vacancy(siteAtomList[0]);
		vacancy.mass=0;
		vacancy.q=0;
		vacancy.dw=0;
		vacancy.Znum=0;
		m_atoms.push_back(vacancy);
		// Probability of a vacancy is the remainder
		probabilities.push_back(1.0-occupancy_sum);
	}
	
	m_alias = alias_method(probabilities);
}

atom &AtomChooser::Choose()
{
	return m_atoms[m_alias.next()];
}

RealVector AtomChooser::GetOccupancyVector()
{
	RealVector occupancies;
	std::vector<atom>::iterator _atom=m_atoms.begin(), end=m_atoms.end();
	for (_atom; _atom!=end; ++_atom)
	{
		occupancies.push_back((*_atom).occ);
	}
	return occupancies;
}

float_tt AtomChooser::SumOccupancy()
{
	float_tt sum=0;
	RealVector probabilities = GetOccupancyVector();
	RealVector::iterator prob=probabilities.begin(), end=probabilities.end();
	for (prob; prob!=end; ++prob) sum+=(*prob);
	return sum;
}