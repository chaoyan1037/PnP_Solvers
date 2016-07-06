// Copyright (c) 2014, ETH Zurich, Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of ETH Zurich or The Regents or University of California
// nor the names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// T. Sattler, C. Sweeney, M. Pollefeys. On Sampling Focal Length Values to
// Solve the Absolute Pose Problem. ECCV 2014.
// Authors: Torsten Sattler (sattlert@inf.ethz.ch),
//          Chris Sweeney (cmsweeney@cs.ucsb.edu)

#ifndef PNP_PARAMS_H_
#define PNP_PARAMS_H_

#include <Eigen/Core>
#include <limits>
#include <vector>
#include <cstdint>

#include "camera_pose.h"
#include "eigen_helpers.h"

using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;

namespace pnp{ 
// Structure encapsulating the parameters of RANSAC.
struct RansacParameters{
  RansacParameters()
    :minimal_sample_number( 3 ),
    failure_probability( 0.01 ),
    min_inlier_ratio( 0.1 ),
    //max_ransac_iterations(std::numeric_limits<uint64_t>::max()),
    max_ransac_iterations( ULLONG_MAX ),
    squared_inlier_threshold( 10.0 ),
    random_seed( -1 ),
    use_T_1_1_test( true ) {
  }
  // minimal number of samples, p3p 3, epnp 5, dlt 6
  int minimal_sample_number;

  // The failure probability of RANSAC. Set to 0.01 means that RANSAC has a 1%
  // chance of missing the correct pose.
  double failure_probability;
  // The minimal assumed inlier ratio, i.e., it is assumed that the given set
  // of correspondences has an inlier ratio of at least min_inlier_ratio.
  // This is required to limit the number of RANSAC iteratios.
  double min_inlier_ratio;
  // Another way to specify the maximal number of RANSAC iterations. In effect,
  // the maximal number of iterations is set to min(max_ransac_iterations, T),
  // where T is the number of iterations corresponding to min_inlier_ratio.
  // This variable is useful if RANSAC is to be applied iteratively, i.e.,
  // first applying RANSAC with an min_inlier_ratio of x, then with one
  // of x-y and so on, and we want to avoid repeating RANSAC iterations.
  // However, the preferable way to limit the number of RANSAC iterations is
  // to set min_inlier_ratio and leave max_ransac_iterations to its default
  // value.
  // Per default, this variable is set to std::numeric_limits<uint64_t>::max().
  uint64_t max_ransac_iterations;
  // The threshold on the squared reprojection error to be used to distinguish
  // between inliers and outliers. Given in pixels squared.
  double squared_inlier_threshold;
  // Whether to use a random seed (-1) or the given seed (>0).
  int random_seed;
  // Whether to use the T_{d,d}, with d=1, test proposed in
  // Chum, O. and Matas, J.: Randomized RANSAC and T(d,d) test, BMVC 2002.
  // After computing the pose, RANSAC selects one match at random and evaluates
  // all poses. If the point is an inlier to one pose, the corresponding pose
  // is rejected. Notice that if the pose solver returns multiple poses, then
  // at most one pose is correct. If the selected match is correct, then only
  // the correct pose will pass the test. Per default, the test is enabled.
  bool use_T_1_1_test;
};

// the method used to solver pnp problem
enum PnPMethod
{
  P3P = 0,
  P3PF = 1,
  EPNP = 2,
  P6P = 3
};

struct PnPResult {
  PnPResult()
    : num_inliers_( 0 ),
    num_generated_random_samples_( 0 ) {
  }

  // The computed pose. The focal length that was used to estimate the pose can
  // be recovered from the intrinsic calibration of the camera pose.
  CameraPose pose_;

  // The number of inliers for that pose.
  int num_inliers_;

  // The number of random samples generated by the p3pf method.
  // For the depth-first-based approaches, this number corresponds to the number
  // of RANSAC iterations. For the breadth-first method, multiplying this number
  // with the number of tested focal length values gives the corresponding
  // number of RANSAC iterations as this approach uses the same random sample
  // for all focal length parameter values during one iteration of p3pf.
  uint64_t num_generated_random_samples_;

  // The points used to generate the pose.
  std::vector<Vector2d> sample_points_2D_;
  std::vector<Vector3d> sample_points_3D_;
};

// Structure encapsulating the parameters of the p3pf approach.
struct P3PfParameters{
  P3PfParameters() {
    focal_length_values_.clear();
    prior_probabilities_.clear();
    cdf_inlier_ratios_.clear();
  }

  // Parameters specific to ransac.
  RansacParameters ransac_parameters_;

  // The set of focal length values to be used for camera pose estimation.
  // For the breadth-first approach, the ordering of the focal length values
  // can be arbitrary. For the depth-first and iterative depth-first variants
  // of p3pf, it is expected that the focal length values are sorted in
  // descending (a-priori) probability of being correct. Giving the focal length
  // values in any other ordering results in non-optimal run-times.
  // For probabilistic p3pf, the focal_length_values should be given in
  // ascending order.
  std::vector<double> focal_length_values_;

  // The a-priori probability for each focal length value. Has to be in the same
  // order as focal_length_values_. This only has to be set if probabilistic
  // p3pf is used.
  std::vector<double> prior_probabilities_;

  // The cumulative distribution function of the inlier ratios (depending on
  // the matching algorithm used to obtain the 2D-3D matches). The algorithm
  // makes the assumption of equidistant values, i.e., if cdf_inlier_ratios_
  // has 11 entries, probabilistic p3pf will interpret it the way that
  // cdf_inlier_ratios_[0] corresponds to an inlier ratio of 0,
  // cdf_inlier_ratios_[1] to an inlier ratio of 0.1, cdf_inlier_ratios_[2] to
  // 0.2, and cdf_inlier_ratios_[10] to 1.0.
  // Only needs to be specified if probabilistic p3pf is used.
  std::vector<double> cdf_inlier_ratios_;
  
};

/*  general pnp( p3p epnp dlt) parameters*/
struct GpnpParameters
{
  GpnpParameters(){ 
    pnp_method_ = P3PF;
    refine_pose = true;
    K_.setIdentity();
  }
  // whether refine pose use all inliers, only for epnp and dlt
  bool refine_pose;

  PnPMethod pnp_method_;

  // Parameters specific to ransac.
  RansacParameters ransac_parameters_;

  // The internal calibration of the camera. K_ has the form diag(f, f, 1),
  // where f is the focal length of the camera.
  Eigen::Matrix3d K_;

};

} // namespace pnp


#endif  // P3PF_PARAMS_H_
