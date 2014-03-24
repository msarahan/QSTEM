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

// This file is an altered version of the Alias method, posted by
// Author: Jaewon Jung (jj@notsharingmy.info)
// at http://ideone.com/gzFzI
//
// This is a C++ port of Keith's Java implementation here:
// http://www.keithschwarz.com/interesting/code/?dir=alias-method
//
// An implementation of the alias method implemented using Vose's algorithm.
// The alias method allows for efficient sampling of random values from a
// discrete probability distribution (i.e. rolling a loaded die) in O(1) time
// each after O(n) preprocessing time.
//
// For a complete writeup on the alias method, including the intuition and
// important proofs, please see the article "Darts, Dice, and Coins: Smpling
// from a Discrete Distribution" at
//
// http://www.keithschwarz.com/darts-dice-coins/
//

#ifndef ATOM_CHOOSER_H
#define ATOM_CHOOSER_H

#include "stemtypes_fftw3.hpp"
#include "random.hpp"

#include <vector>
#include <algorithm>

namespace QSTEM
{

class alias_method
{
public:
/**
* Constructs a new alias_method to sample from a discrete distribution and
* hand back outcomes based on the probability distribution.
* <p>
* Given as input a list of probabilities corresponding to outcomes 0, 1,
* ..., n - 1, along with the random number generator that should be used
* as the underlying generator, this constructor creates the probability
* and alias tables needed to efficiently sample from this distribution.
*
* @param probabilities The list of probabilities.
*/
	alias_method(){};
	alias_method(const RealVector& probability)
	: probability_(probability)
	{
		// Allocate space for the alias table.
		alias_.resize(probability_.size());
 
		// Compute the average probability and cache it for later use.
		const float_tt average = 1.0 / probability_.size();
 
		// two stacks to act as worklists as we populate the tables
		std::vector<int> small, large;
 
		// Populate the stacks with the input probabilities.
		for(size_t i=0; i<probability_.size(); ++i)
		{
			// If the probability is below the average probability, then we add
			// it to the small list; otherwise we add it to the large list.
			if (probability_[i] >= average)
				large.push_back(i);
			else
				small.push_back(i);
		}
 
		// As a note: in the mathematical specification of the algorithm, we
		// will always exhaust the small list before the big list. However,
		// due to floating point inaccuracies, this is not necessarily true.
		// Consequently, this inner loop (which tries to pair small and large
		// elements) will have to check that both lists aren't empty.
		while(!small.empty() && !large.empty())
		{
			// Get the index of the small and the large probabilities.
			int less = small.back(); small.pop_back();
			int more = large.back(); large.pop_back();
 
			alias_[less] = more;
 
			// Decrease the probability of the larger one by the appropriate
			// amount.
			probability_[more] = (probability_[more] + probability_[less]) - average;
 
			// If the new probability is less than the average, add it into the
			// small list; otherwise add it to the large list.
			if (probability_[more] >= average)
				large.push_back(more);
			else
				small.push_back(more);
		}
 
		// At this point, everything is in one list, which means that the
		// remaining probabilities should all be 1/n. Based on this, set them
		// appropriately. Due to numerical issues, we can't be sure which
		// stack will hold the entries, so we empty both.
		while(!small.empty())
		{
			probability_[small.back()] = average;
			small.pop_back();
		}
		while(!large.empty())
		{
			probability_[large.back()] = average;
			large.pop_back();
		}
 
		// These probabilities have not yet been scaled up to be such that
		// 1/n is given weight 1.0. We do this here instead.
		int	n = static_cast<int>(probability_.size());
		std::transform(probability_.cbegin(), probability_.cend(), probability_.begin(),
						[n](float_tt p){ return p * n; });
	}
 
	/**
	* Samples a value from the underlying distribution.
	*
	* @return A random value sampled from the underlying distribution.
	*/
	int next()
	{
		// Generate a fair die roll to determine which column to inspect.
		int column = random_integer() % probability_.size();
 
		// Generate a biased coin toss to determine which option to pick.
		bool coinToss = ran1() < probability_[column];
 
		// Based on the outcome, return either the column or its alias.
		return coinToss? column : alias_[column];
	}
 
private:
	// The probability and alias tables.
	std::vector<int> alias_;
	RealVector probability_;
};


/**
Provide random choices from among atoms (or vacancies) that share a site.

To use this class, instantiate it with a vector of atoms (these should share a site, i.e., 
	have the same x, y, z values), then call the Choose method when you want a random atom.
*/

class AtomChooser
{
public:
	AtomChooser(std::vector<atom> siteAtomList);
	atom &Choose(); // returns a random choice from siteAtomList
private:
	RealVector GetOccupancyVector();
	float_tt SumOccupancy();

	alias_method m_alias; // alias method object, used to choose an integer index of which atom
	std::vector<atom> m_atoms; // list of atoms for this site
};

}
#endif