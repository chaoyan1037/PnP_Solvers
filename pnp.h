
#ifndef PNP_H_
#define PNP_H_

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <string>
#include <vector>
#include <cstdint>
#include <iostream>

#include "pnp_params.h"
#include "camera_pose.h"
#include "p3pf/src/p3pf.h"
#include "epnp/epnp.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using p3pf::P3Pf;

namespace pnp {


struct PnPParameters
{
  //construct function
  PnPParameters();
  ~PnPParameters();

  void InitPnPMethod( const PnPMethod method );

  // indicate current parameter is for which kind of pnp algorithm
  PnPMethod pnp_method;

  // parameter for p3pf
  P3PfParameters* ptr_p3pf_params;
  GpnpParameters* ptr_pnp_params;

};

/* 
  pnp solver consists of p3p(3), p3pf(3), epnp(5), dlt(6), the number in
  the parenthesis indicates the minimal points to solver the problem.
*/
class PnP;
class PnPSolver
{
private:

  PnPMethod pnp_method;
  
  p3pf::P3Pf * ptr_P3Pf;

  // p3p(3), epnp(5), dlt(6)
  PnP * ptr_pnp;

public:

  PnPSolver();
  ~PnPSolver();

  void GetMethod( PnPMethod& method ) const;

  // set the pnp method used to solve the pnp problem
  void Init( const PnPParameters& Parameters );

  void ComputePose(
    const std::vector<Vector2d> &points2D,
    const std::vector<Vector3d> &points3D,
    PnPResult &result ) const;

};
/*
  pnp solver implementation with Ransac loop,  p3p(3), epnp(5), dlt(6)
*/
class PnP{
private:

  GpnpParameters pnp_params_;

  // Do T11Test if T11Test is set
  bool PassesT11Test( 
    const Vector2d& point2D,
    const Vector3d& point3D,
    const CameraPose& pose) const;
  
  bool PassesT11Test(
    const Vector2d& point2D,
    const Vector3d& point3D,
    const std::vector<CameraPose>& vecPose ) const;

  double ComputeSquaredReprojectionError(
    const Vector2d &p_i,
    const Vector3d &p_w,
    const CameraPose &pose) const;

  // find out the number of inliers 
  int EvaluatePose(
    const std::vector<Vector2d> &points2D,
    const std::vector<Vector3d> &points3D,
    const CameraPose &pose,
    std::vector<bool>& inliers ) const;

public:
  PnP(){};
  ~PnP(){};

  // set the pnp method used to solve the pnp problem
  void Init( const GpnpParameters& Parameters );

  void ComputePoseP3P(    
    const std::vector<Vector2d> &points2D,
    const std::vector<Vector3d> &points3D,
    PnPResult &result) const;

  void ComputePoseEpnp(
    const std::vector<Vector2d> &points2D,
    const std::vector<Vector3d> &points3D,
    PnPResult &result ) const;

  void ComputePoseP6P(
    const std::vector<Vector2d> &points2D,
    const std::vector<Vector3d> &points3D,
    PnPResult &result) const;
};

}

#endif