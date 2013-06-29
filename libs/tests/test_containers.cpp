#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "data_containers.h"
#include <iostream>
#include <cstdio>

std::string wavefile = "tests/data/mulswav_16_2.img";
std::string dpfile = "tests/data/diffAvg_0_16.img";
std::string detfile = "tests/data/detector1_2.img";

struct WaveFixture {
  WaveFixture():
    wave(WavePtr( new WAVEFUNC(400, 400, 0.0625, 0.0625)))
  { 
    std::cout << "setup wave fixture" << std::endl; 
  }
  ~WaveFixture()
  { std::cout << "teardown wave fixture" << std::endl; }

  WavePtr wave;
};

struct DetectorFixture {
  DetectorFixture() :
    det(DetectorPtr(new Detector(20, 20, 0.402f, 0.402f)))
  { 
    std::cout << "setup detector fixture" << std::endl; 
  }
  ~DetectorFixture()
  { std::cout << "teardown detector" << std::endl; }

  DetectorPtr det;
};

BOOST_FIXTURE_TEST_SUITE (TestWave, WaveFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  // Check array allocation
  BOOST_CHECK(wave->diffpat != NULL);
  BOOST_CHECK(wave->avgArray != NULL);
  BOOST_CHECK(wave->wave != NULL);
}

// Test image reading (uses a known-good image)
BOOST_AUTO_TEST_CASE (testWaveRead)
{
  wave->ReadWave(wavefile.c_str());
  // test reading the thickness (accuracy within 0.001%.)
  BOOST_CHECK_CLOSE(wave->GetThickness(), 78.0999, 0.001);

  // These tests are of questionable value, because the size and 
  //   resolution should be determined when the wave is created, and should not change from loading a file.
  // test reading the wave array size
  /*
  BOOST_CHECK(wave->GetSizeX()==400);
  BOOST_CHECK(wave->GetSizeY()==400);
  BOOST_CHECK_CLOSE(wave->GetResolutionX(), 0.0625, 0.001);
  BOOST_CHECK_CLOSE(wave->GetResolutionY(), 0.0625, 0.001);
  */
  
}


// Test image saving
// saves the data from the known good image, then reads it back in.  
//    This test cannot pass if the image reading test does not pass.
BOOST_AUTO_TEST_CASE (testWaveSave)
{
  wave->ReadWave(wavefile.c_str());
  wave->WriteWave("wave_output_test.img");
  wave->ReadWave("wave_output_test.img");
  BOOST_CHECK_CLOSE(wave->GetThickness(), 78.0999, 0.001);
  remove( "wave_output_test.img" );
}

BOOST_AUTO_TEST_SUITE_END( )


BOOST_FIXTURE_TEST_SUITE (TestDetector, DetectorFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  // Check array allocation
  BOOST_CHECK(det->image != NULL);
  BOOST_CHECK(det->image2 != NULL);
}

BOOST_AUTO_TEST_CASE (testSetComment)
{
	
}

BOOST_AUTO_TEST_SUITE_END( )
