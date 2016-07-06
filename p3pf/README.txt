P3P(f) on Windows VS2013.

The original software is implemented on Linux, and have been tested on Mac OS, 
however it need some modifications to work on Windows, both source code and 
CMakeLists.txt.
Specifically, I modify the CMakeLists.txt so that CMake automatically include
and link the glog and gflags libraries. For details, please refer to the CMakeLists.

For your usage, you need only modify the Eigen include path.
--------------------------------------------------------------------------------



P3P(f) Version 1.0
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
What is this?
--------------------------------------------------------------------------------
This is the reference implementation of the P3P(f) method that estimates the
pose of a camera with unknown focal length. Please see the original publication
  Torsten Sattler, Chris Sweeney, Marc Pollefeys.
  On Sampling Focal Length Values to Solve the Absolute Pose Problem.
  ECCV 2014.
for a detailed description of the method.

--------------------------------------------------------------------------------
License and conditions of usage
--------------------------------------------------------------------------------
The software in the src/ directory is licensed under the BSD 3-Clause License
(see http://opensource.org/licenses/BSD-3-Clause and the LICENSE file).

The software package also contains a version of Laurent Kneip's code to estimate
the camera pose of a calibrated camera from three 2D-3D matches, located in
ext/p3p_code_kneip/. The original software has been modified to use the Eigen
library and is licensed under the BSD-new license. Please see
ext/p3p_code_kneip/Readme.txt for a description of the method and the conditions
of usage. These files have been provided for convenience only.

The software package also contains cmake files taken from the Ceres library. The
related copyright notices and licenses are contained in the corresponding cmake
files. These files have been provided for convenience only.

The software package also contains some files released under the BSD-new license
that are copyright of The Regents of the University of California.

If you are using the P3P(f) implementation for research, we kindly ask you to
cite the paper
  Torsten Sattler, Chris Sweeney, Marc Pollefeys.
  On Sampling Focal Length Values to Solve the Absolute Pose Problem.
  ECCV 2014.
in your publication.

--------------------------------------------------------------------------------
Installation
--------------------------------------------------------------------------------
Requirements:
* CMake (http://www.cmake.org/)
* Eigen 3.2.0 or newer (http://eigen.tuxfamily.org/)
* glog (http://code.google.com/p/google-glog/)
* C+11

We have tested this software under Ubuntu 14.04 and Mac OS X Mavericks.
We only support these two systems, but compiling the software under Windows
should be simple as well.

To build and install the software, go to the directory containing the README.txt
file and type

mkdir build
mkdir release
cd build
cmake -DCMAKE_INSTALL_PREFIX=../release ..
make install

This will install the P3P(f) library under release/lib, the header files under
release/include, and an executable for an example program (p3pf_example) under
release/bin.

To use the example program, compile the code using the steps above and navigate
to the bin directory. Run the command

./p3pf_example --logtostderr

and this will run p3pf on a synthetic test example and log the results.

Per default, the software is build in Release mode.

--------------------------------------------------------------------------------
Using the P3P(f) software
--------------------------------------------------------------------------------
Please see src/p3pf_example.cc for an example on how to use the P3P(f) algorithm.
In this file, random 2D-3D matches are generated and the pose is estimated
using the P3P(f) method. The file also contains the priors on the opening angles
used in the paper.

For further questions on using the software, please contact Torsten Sattler
(torsten.sattler@inf.ethz.ch).

Please understand that we do not have the resources to provide full technical
support.
