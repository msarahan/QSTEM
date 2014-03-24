#ifndef IMAGELIB_H
#define IMAGELIB_H

#include "stemtypes_fftw3.hpp"
#include <vector>
#include <map>
#include <string>
#include "boost/shared_ptr.hpp"
#include "data_IO/data_io_factories.hpp"
#include "crystal.hpp"

/**************************************************************
 * Here is how to use the new image writing routines
 *
 * ImageIOPtr imageio = ImageIOPtr(new CImageIO(nx,ny))
 *
 * imageIO->WriteImage((void**)real_image,filename);
 * imageIO->WriteImage((void**)complex_image,filename);
 *
 **************************************************************
 * Reading an image works like this:
 * 
 * ImageIOPtr imageio = ImageIOPtr(new CImageIO(nx,ny))
 * imageio->ReadImage((void **)pix,nx,ny,fileName);
 *
 * Note that header parameters are persistent on any object.  You
 *   should use the various Set* functions to set parameters as
 *   necessary.  You should not need to read values from this class - 
 *   only set them.  They will be recorded to any file saved from this
 *   this object.
 **************************************************************/

namespace QSTEM
{

/**
Helper class for image input/output.  This class stores an image size and initializes 
reader/writer classes appropriately.  It provides several API options depending on what exact details
you want to pass to the reader/writer plugin.
*/
class CImageIO {
  int m_nx,m_ny;
  char m_buf[200];  // General purpose temporary text buffer
  DataWriterPtr m_imageWriter;
  DataReaderPtr m_imageReader;
public:
  CImageIO(int nx, int ny, const std::string input_extension=".img", 
           const std::string output_extension=".img");

  void CreateRealDataSet(const std::string &name, const std::vector<unsigned int> &positions);
  void CreateComplexDataSet(const std::string &name, const std::vector<unsigned int> &positions);
  
  void WriteImage(const RealVector &pix, const std::string &fileName, 
                      std::map<std::string, double> &params,
                      const std::string &comment, 
                      const std::vector<unsigned> &position);
  // If you want a comment and parameters, but no position
  inline void WriteImage(const RealVector &pix, const std::string &fileName, std::map<std::string, double> &parameters, 
                             const std::string &comment)
  {
    std::vector<unsigned> position;
    WriteImage(pix, fileName, parameters, comment, position);
  }
  // If you want to add a comment, but no parameters
  inline void WriteImage(const RealVector &pix, const std::string &fileName, const std::string &comment,
                             const std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    WriteImage(pix, fileName, params, comment, position);
  }
  // If you want to add parameters, but no comment
  inline void WriteImage(const RealVector &pix, const std::string &fileName, std::map<std::string, double> &params,
                               const std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::string comment=std::string();
    WriteImage(pix, fileName, params, comment, position);
  }  
  // If you don't care about parameters or comment, use this simplified overload:
  inline void WriteImage(const RealVector &pix, const std::string &fileName, 
                        const std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    std::string comment=std::string();
      
    WriteImage(pix, fileName, params, comment, position);
  }


  void WriteImage(const ComplexVector &pix, const std::string &fileName,
                         std::map<std::string, double> &params,
                         const std::string &comment,
                         const std::vector<unsigned> &position);
  // If you want a comment and parameters, but no position
  inline void WriteImage(const ComplexVector &pix, const std::string &fileName, std::map<std::string, double> &parameters, 
                             const std::string &comment)
  {
    std::vector<unsigned> position;
    WriteImage(pix, fileName, parameters, comment, position);
  }

  // If you want to add a comment, but no parameters
  inline void WriteImage(const ComplexVector &pix, const std::string &fileName, const std::string &comment,
                             const std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    WriteImage(pix, fileName, params, comment, position);
  }
  // If you want to add parameters, but no comment
  inline void WriteImage(const ComplexVector &pix, const std::string &fileName, std::map<std::string, double> &params,
                               const std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::string comment=std::string();

    WriteImage(pix, fileName, params, comment, position);
  }  
  // If you don't care about parameters or comment, use this simplified overload:
  inline void WriteImage(const ComplexVector &pix, const std::string &fileName, 
                        const std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    std::string comment=std::string();

    WriteImage(pix, fileName, params, comment, position);
  }


  void ReadImage(RealVector &pix, std::string &fileName, std::map<std::string, double> &params,
                 std::string &comment,
                 std::vector<unsigned> position=std::vector<unsigned>());
  void ReadImage(ComplexVector &pix, std::string &fileName, std::map<std::string, double> &params,
                 std::string &comment,
                 std::vector<unsigned> position=std::vector<unsigned>());
  // If you don't care about parameters, this reads the data without you passing them in.
  template <typename T>
  inline void ReadImage(T &pix, std::string &fileName, std::string &comment,
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    ReadImage(pix, fileName, params, comment, position);
  }
  // If you want the parameters, but don't care about the comment
  template <typename T>
  inline void ReadImage(T &pix, std::string &fileName, std::map<std::string, double> &params,
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::string comment;
    ReadImage(pix, fileName, params, comment, position);
  }
  // If you don't care about comment or parameters
  template <typename T>
  inline void ReadImage(T &pix, std::string &fileName,
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    std::string comment;
    ReadImage(pix, fileName, params, comment, position);
  }
  
private:
  std::vector<unsigned> GetShapeVector();
};

typedef boost::shared_ptr<CImageIO> ImageIOPtr;

}  // end namespace QSTEM

#endif
