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

#ifndef P3PF_H_
#define P3PF_H_

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <string>
#include <vector>

#include "eigen_helpers.h"
#include "pnp_params.h"
//#include "src/p3pf_params.h"

namespace p3pf {
using Eigen::Vector2d;
using Eigen::Vector3d;
using pnp::CameraPose;
using pnp::PnPResult;
using pnp::P3PfParameters;

class P3PSolver;


// Implements a method that estimates the camera pose by combining the use of a
// 3 point pose solver inside a normal RANSAC loop with sampling focal length
// values from a discrete set of parameters.
// In each RANSAC iteration, this class randomly selects a focal length
// depending on a probability distribution that is updated as the algorithm
//  progresses.
class P3Pf {
 public:
  P3Pf();
  virtual ~P3Pf();

  // Initializes the P3PfBreadthFirst instance using the given parameters and
  // also initializes the three point pose solver. CHECKs that
  // parameters.focal_length_values_ and parameters.prior_probabilities_ are not
  // empty and have the same length. CHECKs that parameters.cdf_inlier_ratios_
  // is not empty.
  virtual void Init(const P3PfParameters &parameters);

  // Estimates the camera pose from a given set of 2D-3D correspondences by
  // combining a 3 point pose solver inside a RANSAC loop with focal length
  // parameter value sampling. In each step, one of the focal length values
  // from parameters.focal_length_values is selected randomly according to
  // some probability distribution that is updated as the algorithm
  // progresses. The best pose, i.e., the pose with the largest number of 3D
  // points is returned in the result parameter.
  virtual void ComputePose(
      const std::vector<Vector2d> &points2D,
      const std::vector<Vector3d> &points3D,
      PnPResult *result) const;

 protected:
  // Evaluates a single point per the T{1,1} test and returns true if that point
  // is an inlier and false otherwise. This is used to speed up the evaluation
  // of correspondences.
  /*
    bool PassesT11Test(const Vector3d &points3d,
                     const Vector2d &points2d,
                     P3PSolver *solver) const;
  */

  // Updates the sampling probabilities by calling the independent and dependent
  // methods below as necessary.
  double UpdateSamplingProbabilities(
      const std::vector<uint64_t> &num_samples_per_value,
      const int best_focal_length_value,
      const double best_inlier_ratio,
      std::vector<double> *probabilities,
      bool *better_inlier_ratio_possible) const;

  // Given the number of samples taken for every focal length value and the
  // best inlier ratio found so far, computes the probability for each focal
  // length with which it will be selected in the next iteration. Returns a
  // vector containing the probability for each focal length. Assumes that the
  // individual focal length values are independent from each other.
  // For each focal length value, computes the probability of finding a better
  // model according to the cdf of the inlier ratios and returns the maximum
  // among these probabilities.
  double ComputeSamplingProbabilitiesIndependent(
      const std::vector<uint64_t> &num_samples_per_value,
      const double best_inlier_ratio,
      std::vector<double> *probabilities,
      bool *better_inlier_ratio_possible) const;

  // Given the number of samples taken for every focal length value, the best
  // inlier ratio found so far, and the index of the focal length for which the
  // best model was found, computes the probability for each focal
  // length with which it will be selected in the next iteration. Returns a
  // vector containing the probability for each focal length. Assumes that the
  // is a dependence between the individual focal length values.
  // CHECKs that the index of the best focal length value is in the range
  // [ 0, parameters.focal_length_values_.size() ).
  // For each focal length value, computes the probability of finding a better
  // model according to the cdf of the inlier ratios and returns the maximum
  // among these probabilities.
  double ComputeSamplingProbabilitiesDependent(
      const std::vector<uint64_t> &num_samples_per_value,
      const int best_focal_length_value,
      const double best_inlier_ratio,
      std::vector<double> *probabilities,
      bool *better_inlier_ratio_possible) const;

  // Given the probabilities of the focal length values and a number chosen
  // uniformly at random from the interval [0, 1), returns the index of the
  // focal length corresponding to this value according to the probability
  // distribution defined by probabilities. Assumes that probabilities has the
  // same length as parameters.focal_length_values_.
  int SampleFocalLength(const std::vector<double> &probabilities,
                        const double random_value) const;

  // Given a number of iterations k, returns the inlier ratio such that
  // the probability of having missed an all-inlier sample with in k steps is
  // at most parameters.ransac_parameters.failure_probability.
  double GetInlierRatioFromNumberIterations(
      const uint64_t num_iterations) const;

  // The multiplicative factor used to discretize the inlier ratio so that we
  // can perform a look-up in parameters.cdf_inlier_ratios_.
  // Is initialized in Init().
  double scaling_factor_inlier_ratios_;

  // Look-up table that stores for each number of iterations the corresponding
  // inlier ratio depending on parameters.ransac_parameters_. Initialized
  // in Init(). Notice that this might lead to large memory consumptions in
  // case that parameters.ransac_parameters_.min_inlier_ratio is extremely low.
  std::vector<double> inlier_ratio_per_num_iterations_;

  P3PfParameters params_;
};

}  // namespace p3pf

#endif  // P3PF_H_
