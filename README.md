# PnP_Solvers
various PnP slover implementations, include p3p, p3pf, epnp, dlt

the pnp_solvers is based on the p3pf originally provided by Torsten Sattler
(http://people.inf.ethz.ch/sattlert/). I add another two pnp algorithms, epnp 
and dlt.

Dependency:
OpenCV(epnp), Eigen, and Google gflags and glogs(For this part you can refer
to the README.txt in ./p3pf).

Notes:

CMakeLists.txt is provided and you can complie the program with CMake.
You need modify the CMakeLists.txt to find OpenCV library, also the Eigen include path. 
pnp_example.cc file are provided to demonstrate how to use the PnP_Solvers. 



Reference Papers and Books: 

Torsten Sattler, Chris Sweeney, Marc Pollefeys.
On Sampling Focal Length Values to Solve the Absolute Pose Problem.
ECCV 2014.

Laurent Kneip, Davide Scaramuzza, Roland Siegwart 
A Novel Parametrization of the Perspective-Three-Point Problem for 
a Direct Computation of Absolute Camera Position and Orientation.
CVPR 2011.

Lepetit Vincent, Francesc Morenonoguer, Pascal Fua
EPnP: An Accurate O(n) Solution to the PnP Problem.
IJCV 2009

Hartley Richard, Andrew Zisserman
Multiple View Geometry
IJCV 2000