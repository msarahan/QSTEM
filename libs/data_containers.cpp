#include "stdio.h"
#include <string.h>
#include "data_containers.h"

WAVEFUNC::WAVEFUNC(int x, int y) :
detPosX(0),
detPosY(0),
iPosX(0),
iPosY(0),
thickness(0.0),
nx(0),
ny(0)
{

  std::string waveFile;
  std::string waveFileBase = "mulswav";
  nx = x;
  ny = y;
  diffpat = QSfMat::Zero(nx, ny);
  avgArray = QSfMat::Zero(nx, ny);
  wave = QScMat::Zero(nx, ny);
  
  /*
#if FLOAT_PRECISION == 1
	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
		fftMeasureFlag);
	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
		fftMeasureFlag);
#endif
  */
        waveFile = waveFileBase+".img";
        fileout = waveFile;
        fileStart = "mulswav.img";
}

WAVEFUNC::WAVEFUNC( WAVEFUNC& other )
{
	iPosX = other.iPosX;
	iPosY = other.iPosY;
	detPosX = other.detPosX;
	detPosY = other.detPosY;

        fileStart = other.fileStart;
        fileout = other.fileout;
        avgName = other.avgName;

	nx = other.nx;
	ny = other.ny;

	diffpat = QSfMat::Zero(nx, ny);
	avgArray = QSfMat::Zero(nx, ny);

    wave = QScMat::Zero(nx,ny);

#if FLOAT_PRECISION == 1
        //	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
        //	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
        //	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
        //		fftMeasureFlag);
        //	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
        //		fftMeasureFlag);
#endif
}

// reset the wave's thickness
void WAVEFUNC::ZeroWave(void)
{

}

MULS::MULS () 
{
	initMuls();
};

MULS::MULS (int slices)  
{
	// muls.areaAIS = 1.0;

	/* make multislice read the inout files and assign transr and transi: */
	//muls.trans = NULL;
	//muls.cz = NULL;  // (float_tt *)malloc(muls.slices*sizeof(float_tt));

	//muls.kx = NULL;
	//muls.kx2= NULL;
	//muls.ky = NULL;
	//muls.ky2= NULL;

	/****************************************************/
	/* copied from slicecell.c                          */
	//muls.pendelloesung = NULL;
}

void MULS::initMuls(int slices)
{
	initMuls();
	int sCount;
	/* general setup: */
	// cin2 should be a vector of strings - with each element being "a" followed by the slice index.  
	//     Is a a character (integer), and this is some kind of offset based on that?
	for (sCount =0;sCount<slices;sCount++)
		cin2[sCount] = 'a'+sCount;
	for (sCount = slices;sCount < NCINMAX;sCount++)
		cin2[sCount] = 0;

}

void MULS::initMuls()
{
	atomRadius= 5.0; /* radius in A for making the potential boxes */
	saveFlag=0;
	sigmaf=0; dfdelt=0; acmax=0; acmin=0; aobj=0; Cs=0; aAIS=0;
	// Tomography stuff - tomoCount = 0 indicates no Tomo simulation.
	tomoTilt=0; tomoStart=0; tomoStep=0; tomoCount=0;
	onlyFresnel=0;
	czOffset=0; /* defines the offset for the first slice in 
						fractional coordinates        */
	normHolog=0;
	gaussianProp=0;

	sparam = QSfVec::Zero(NPARAM);
}

ElTable::ElTable()
{
	elements = std::vector<std::string>(104);
	elements.push_back("H");
	elements.push_back("He");
	elements.push_back("Li");
	elements.push_back("Be");
	elements.push_back("B");
	elements.push_back("N");
	elements.push_back("C");
	elements.push_back("O");
	elements.push_back("F");
	elements.push_back("Ne");
	elements.push_back("Na");
	elements.push_back("Mg");
	elements.push_back("Al");
	elements.push_back("Si");
	elements.push_back("P");
	elements.push_back("S");
	elements.push_back("Cl");
	elements.push_back("Ar");
	elements.push_back("K");
	elements.push_back("Ar");
	elements.push_back("K");
	elements.push_back("Ca");
	elements.push_back("Sc");
	elements.push_back("Ti");
	elements.push_back("V");
	elements.push_back("Cr");
	elements.push_back("Mn");
	elements.push_back("Fe");
	elements.push_back("Co");
	elements.push_back("Ni");
	elements.push_back("Cu");
	elements.push_back("Zn");
	elements.push_back("Ga");
	elements.push_back("Ge");
	elements.push_back("As");
	elements.push_back("Se");
	elements.push_back("Br");
	elements.push_back("Kr");
	elements.push_back("Rb");
	elements.push_back("Sr");
	elements.push_back("Y");
	elements.push_back("Zr");
	elements.push_back("Nb");
	elements.push_back("Mo");
	elements.push_back("Tc");
	elements.push_back("Ru");
	elements.push_back("Rh");
	elements.push_back("Pd");
	elements.push_back("Ag");
	elements.push_back("Cd");
	elements.push_back("In");
	elements.push_back("Sn");
	elements.push_back("Sb");
	elements.push_back("Te");
	elements.push_back("I");
	elements.push_back("Xe");
	elements.push_back("Cs");
	elements.push_back("Ba");
	elements.push_back("La");
	elements.push_back("Ce");
	elements.push_back("Pr");
	elements.push_back("Nd");
	elements.push_back("Pm");
	elements.push_back("Sm");
	elements.push_back("Eu");
	elements.push_back("Gd");
	elements.push_back("Tb");
	elements.push_back("Dy");
	elements.push_back("Ho");
	elements.push_back("Er");
	elements.push_back("Tm");
	elements.push_back("Yb");
	elements.push_back("Lu");
	elements.push_back("Hf");
	elements.push_back("Ta");
	elements.push_back("W");
	elements.push_back("Re");
	elements.push_back("Os");
	elements.push_back("Ir");
	elements.push_back("Pt");
	elements.push_back("Au");
	elements.push_back("Hg");
	elements.push_back("Tl");
	elements.push_back("Pb");
	elements.push_back("Bi");
	elements.push_back("Po");
	elements.push_back("At");
	elements.push_back("Rn");
	elements.push_back("Fr");
	elements.push_back("Ra");
	elements.push_back("Ac");
	elements.push_back("Th");
	elements.push_back("Pa");
	elements.push_back("U");
	elements.push_back("Np");
	elements.push_back("Pu");
	elements.push_back("Am");
	elements.push_back("Cm");
	elements.push_back("Bk");
	elements.push_back("Cf");
	elements.push_back("Es");
	elements.push_back("Fm");
	elements.push_back("Md");
	elements.push_back("No");
	elements.push_back("Lr");
}