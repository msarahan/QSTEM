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

#ifndef ELTABLE_H
#define ELTABLE_H

#include <string>
#include <boost/algorithm/string.hpp>

static std::string elTable = 
	"H HeLiBeB C N O F NeNaMgAlSiP S Cl"
	"ArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBr"
	"KrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTe"
	"I XeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
	"YbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"
	"FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLr";

inline int getZNumber(std::string element) 
{  
	boost::algorithm::trim(element);
  if (element.size()==1)
    {
      element += " ";
    }
  size_t elem = elTable.find(element);
  if (elem == std::string::npos)
    return 0;	
  else
    return (int)(elem)/2+1;

}

inline std::string getSymbol(unsigned ZNumber)
{
  if (ZNumber==0) return "";
  std::string elem = elTable.substr(2*ZNumber-2, 2);
  if (elem[1] == ' ') elem.erase(1,1);
  return elem;
}

#endif
