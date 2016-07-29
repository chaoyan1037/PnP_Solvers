# PnP_Solvers
Various PnP slover implementations, include p3p, p3pf, epnp, dlt

the pnp_solvers is based on the p3pf originally provided by [Torsten Sattler](http://people.inf.ethz.ch/sattlert/). 
I add another two pnp algorithms, epnp and dlt. And they work in the same framework, so you can easily change the method you want to use.

## Dependency:
OpenCV(epnp), Eigen, and Google gflags and glogs(For this part you can refer
to the README.txt in ./p3pf).

## Notes:

CMakeLists.txt is provided and you can complie the program with CMake.
You need modify the CMakeLists.txt to find OpenCV library, also the Eigen include path. 
Serval pnp_example.cc files are provided to demonstrate how to use the PnP_Solvers. 



## Reference Papers and Books: 

Torsten Sattler, Chris Sweeney, Marc Pollefeys.<br> 
On Sampling Focal Length Values to Solve the Absolute Pose Problem.<br> 
ECCV 2014.

Laurent Kneip, Davide Scaramuzza, Roland Siegwart <br> 
A Novel Parametrization of the Perspective-Three-Point Problem for 
a Direct Computation of Absolute Camera Position and Orientation.<br> 
CVPR 2011.<br> 

Lepetit Vincent, Francesc Morenonoguer, Pascal Fua<br> 
EPnP: An Accurate O(n) Solution to the PnP Problem.<br> 
IJCV 2009<br> 

Hartley Richard, Andrew Zisserman<br> 
Multi View Geometry in Computer Vision.<br> 
