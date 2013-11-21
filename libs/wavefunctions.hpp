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

#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "stemtypes_fftw3.hpp"
#include "imagelib_fftw3.hpp"
#include "memory_fftw3.hpp"
#include "config_readers.hpp"
#include "potential.hpp"

void CreateWaveFunctionDataSets(unsigned x, unsigned y, std::vector<unsigned> positions, std::string output_ext);

static std::string waveFilePrefix="mulswav";
static std::string dpFilePrefix="diff";
static std::string avgFilePrefix="diffAvg";
static std::string probeFilePrefix="probe_wave";
static std::string imageFilePrefix="image";
static std::string waveIntensityFilePrefix="waveIntensity";

// a structure for a probe/parallel beam wavefunction.
// Separate from mulsliceStruct for parallelization.
class WAVEFUNC 
{
  // shared pointer to 
  ImageIOPtr m_imageIO;
public:
  std::string m_fileStart;
  std::string m_avgName;
  std::string m_fileout;
  unsigned m_detPosX, m_detPosY;
  unsigned m_iPosX,m_iPosY;      /* integer position of probe position array */
  unsigned m_nx, m_ny;			/* size of diffpat arrays */
  float_tt **m_diffpat;
  float_tt **m_avgArray;
  float_tt m_thickness;
  float_tt m_intIntensity;
  float_tt m_electronScale;
  float_tt m_beamCurrent;
  float_tt m_dwellTime;
  float_tt m_v0;
  std::vector<unsigned> m_position;
  std::map<std::string, double> m_params;

  std::vector<float_tt> m_kx2,m_ky2,m_kx,m_ky;
  float_tt m_k2max;

  // These are not used for anything aside from when saving files.
  float_tt m_resolutionX, m_resolutionY;

  complex_tt  **m_wave; /* complex wave function */

#if FLOAT_PRECISION == 1
  fftwf_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#else
  fftw_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#endif

public:
  // initializing constructor:
  WAVEFUNC(unsigned nx, unsigned ny, float_tt resX, float_tt resY, std::string input_ext, std::string output_ext);
  WAVEFUNC(ConfigReaderPtr &configReader);
  // define a copy constructor to create new arrays
  WAVEFUNC( WAVEFUNC& other );

  void CreateDataSets();
  void Transmit(PotPtr pot, unsigned sliceIdx);
  void Propagate();
  void FormProbe();

  inline void WriteProbe()
  {
    _WriteWave(probeFilePrefix);
  }
  // WriteImage is for TEM mode
  inline void WriteImage()
  {
    _WriteWave(imageFilePrefix, "Image intensity");
  }
  inline void WriteWaveIntensity()
  {
    _WriteDiffPat(waveIntensityFilePrefix, "Wave intensity");
  }

  inline void WriteWave(std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    m_position.clear();
    _WriteWave(waveFilePrefix, comment, params);
  }
  inline void WriteWave(unsigned navg, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(navg);
    _WriteWave(waveFilePrefix, comment, params);
  }
  inline void WriteWave(unsigned posX, unsigned posY, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(posX, posY);
    _WriteWave(waveFilePrefix, comment, params);
  }

  inline void WriteDiffPat(std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    m_position.clear();
    _WriteDiffPat(dpFilePrefix, comment, params);
  }
  inline void WriteDiffPat(unsigned navg, std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(navg);
    _WriteDiffPat(dpFilePrefix, comment, params);
  }
  inline void WriteDiffPat(unsigned posX, unsigned posY, std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(posX, posY);
    _WriteDiffPat(dpFilePrefix, comment, params);
  }

  inline void WriteAvgArray(std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    m_position.clear();
    _WriteAvgArray(avgFilePrefix, comment, params);
  }
  inline void WriteAvgArray(unsigned navg, std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(navg);
    _WriteAvgArray(avgFilePrefix, comment, params);
  }
  inline void WriteAvgArray(unsigned posX, unsigned posY, std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(posX, posY);
    _WriteAvgArray(avgFilePrefix, comment, params);
  }

  // ReadImage is for TEM mode
  void ReadImage();
  void ReadWave();
  void ReadWave(unsigned navg);
  void ReadWave(unsigned posX, unsigned posY);
  void ReadDiffPat();
  void ReadDiffPat(unsigned navg);
  void ReadDiffPat(unsigned posX, unsigned posY);
  void ReadAvgArray();
  void ReadAvgArray(unsigned navg);
  void ReadAvgArray(unsigned posX, unsigned posY);

private:
  void _WriteWave(std::string &prefix, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>());
  void _WriteDiffPat(std::string &prefix, std::string comment="Diffraction Pattern",
                    std::map<std::string, double>params = std::map<std::string, double>());
  void _WriteAvgArray(std::string &prefix, std::string comment="Average Array",
                     std::map<std::string, double>params = std::map<std::string, double>());

  // For CBED ( &TEM? )
  void SetWavePosition(unsigned navg);
  // For STEM
  void SetWavePosition(unsigned posX, unsigned posY);

};

typedef boost::shared_ptr<WAVEFUNC> WavePtr;

#endif
