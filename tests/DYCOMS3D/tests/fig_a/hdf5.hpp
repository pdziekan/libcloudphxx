#pragma once

#include <blitz/array.h>
#include <H5Cpp.h>
#include <map>

std::map<std::string, double> h5n(
  const string &file
)
{
  H5::H5File h5f(file + "/const.h5", H5F_ACC_RDONLY);
  hsize_t n[3];
  std::map<std::string, double> map;

  h5f.openDataSet("T").getSpace().getSimpleExtentDims(n, NULL);
  map["t"] = n[0];

  {
    auto h5d = h5f.openDataSet("T");

    float dt;
    {
      auto attr = h5d.openAttribute("dt");
      attr.read(attr.getDataType(), &dt);
    }
    map["dt"] = dt;

    auto h5s = h5d.getSpace();
    const hsize_t two = 2, zero = 0;
    float tmp[2];
    h5s.selectHyperslab( H5S_SELECT_SET, &two, &zero);
    h5d.read(tmp, H5::PredType::NATIVE_FLOAT, H5::DataSpace(1, &two), h5s);
    map["outfreq"] = (tmp[1] - tmp[0]) / dt;
  }

  h5f.openDataSet("X").getSpace().getSimpleExtentDims(n, NULL); // X gives cell-border coordinates (+1)
  map["x"] = n[0]-1;
  map["y"] = n[1]-1;
  map["z"] = n[2]-1;

  // read dx,dy,dz
  H5::DataSet h5d = h5f.openDataSet("X");
  H5::DataSpace h5s = h5d.getSpace();

  if (h5s.getSimpleExtentNdims() != 3) 
    error_macro("need 3 dimensions")

  enum {x,y,z};
  h5s.getSimpleExtentDims(n, NULL);

  blitz::Array<float, 3> tmp(n[x], n[y], n[z]);

  hsize_t 
    cnt[3] = { n[x], n[y], n[z] }, 
    off[3] = { 0,    0,    0    };
  h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);

  hsize_t ext[3] = {
    hsize_t(tmp.extent(0)), 
    hsize_t(tmp.extent(1)), 
    hsize_t(tmp.extent(2))
  };
  h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, ext), h5s);
  map["dx"] = tmp(1,0,0) - tmp(0,0,0);
  h5d = h5f.openDataSet("Y");
  h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, ext), h5s);
  map["dy"] = tmp(0,1,0) - tmp(0,0,0);
  h5d = h5f.openDataSet("Z");
  h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, ext), h5s);
  map["dz"] = tmp(0,0,1) - tmp(0,0,0);

  return map;
}

auto h5load(
  const string &file, 
  const string &dataset,
  int at
) -> decltype(blitz::safeToReturn(blitz::Array<float, 3>() + 0))
 {
  notice_macro("about to open file: " << file)
  H5::H5File h5f(file + "/timestep" + zeropad(at, 10) + ".h5", H5F_ACC_RDONLY);

  notice_macro("about to read dataset: " << dataset)
  H5::DataSet h5d = h5f.openDataSet(dataset);
  H5::DataSpace h5s = h5d.getSpace();

  if (h5s.getSimpleExtentNdims() != 3) 
    error_macro("need 3 dimensions")

  hsize_t n[3];
  enum {x,y,z};
  h5s.getSimpleExtentDims(n, NULL);

  blitz::Array<float, 3> tmp(n[x], n[y], n[z]);

  hsize_t 
    cnt[3] = { n[x], n[y], n[z] }, 
    off[3] = { 0,    0,    0    };
  h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);

  hsize_t ext[3] = {
    hsize_t(tmp.extent(0)), 
    hsize_t(tmp.extent(1)), 
    hsize_t(tmp.extent(2))
  };
  h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, ext), h5s);

  return blitz::safeToReturn(tmp + 0);
}
