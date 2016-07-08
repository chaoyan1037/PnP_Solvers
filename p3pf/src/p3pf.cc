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

#include "src/p3pf.h"

#ifndef GLOG_NO_ABBREVIATED_SEVERITIES
#define GLOG_NO_ABBREVIATED_SEVERITIES
#endif
#if _WIN32
#define NOMINMAX
#include <windows.h>
#include <windef.h>
#endif
#include <glog/logging.h>

#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "camera_pose.h"
#include "src/p3p_solver.h"

using std::max;

namespace p3pf {
using Eigen::Vector2d;
using Eigen::Vector3d;
using pnp::CameraPose;
using pnp::PnPResult;

namespace {

// Quantizes a double value from [0, 1] to the integer range [0, max_val] by
// interpolating the value.
double QuantizeValue(const std::vector<double> &values,
                     const double in_val,
                     const double max_val) 
{
  const double scaled_val = in_val * max_val;
  const int scaled_floor_val = static_cast<int>(std::floor(scaled_val));
  if (scaled_floor_val == values.size() - 1) {
    return values[scaled_floor_val];
  }
  const int scaled_ceil_val = scaled_floor_val + 1;
  const double slope = values[scaled_ceil_val] - values[scaled_floor_val];
  const double x_diff = scaled_val - scaled_floor_val;

  return values[scaled_floor_val] + x_diff * slope;
}

// Returns the maximal number of iterations needed such that the probability
// of having missed the correct match falls below the given threshold.
uint64_t GetMaxRANSACIterations(const double inlier_ratio,
                                const double log_threshold,
                                const bool use_T_1_1_test) 
{
  uint64_t max_iters = 0;
  if (inlier_ratio < 1.0) {
    double prob_all_inliers = inlier_ratio * inlier_ratio * inlier_ratio;
    if (use_T_1_1_test) {
      prob_all_inliers *= inlier_ratio;
    }
    const double num_iters = log_threshold / log(1.0 - prob_all_inliers);
    max_iters = static_cast<uint64_t>(std::floor(num_iters));
  }
  return max_iters + 1;
}

}  // namespace

P3Pf::P3Pf() {}

P3Pf::~P3Pf() {}

void P3Pf::Init(const P3PfParameters &parameters) {
  CHECK_LE(parameters.ransac_parameters_.min_inlier_ratio, 1.0);
  CHECK_GT(parameters.ransac_parameters_.min_inlier_ratio, 0.0);
  CHECK_LT(parameters.ransac_parameters_.failure_probability, 1.0);
  CHECK_GT(parameters.ransac_parameters_.failure_probability, 0.0);
  CHECK(!parameters.focal_length_values_.empty());

  params_.focal_length_values_.assign(parameters.focal_length_values_.begin(),
                                      parameters.focal_length_values_.end());
  params_.prior_probabilities_.assign(parameters.prior_probabilities_.begin(),
                                      parameters.prior_probabilities_.end());
  params_.cdf_inlier_ratios_.assign(parameters.cdf_inlier_ratios_.begin(),
                                    parameters.cdf_inlier_ratios_.end());
  params_.ransac_parameters_ = parameters.ransac_parameters_;

  CHECK_EQ(parameters.focal_length_values_.size(),
           parameters.prior_probabilities_.size());
  CHECK_GE(parameters.cdf_inlier_ratios_.size(), 2);

  // Computes the scaling factor that transforms an inlier ratio value from
  // [0, 1] to a position in the look-up table for the cdf.
  scaling_factor_inlier_ratios_ = static_cast<double>(parameters
      .cdf_inlier_ratios_.size() - 1);

  // Initializes the look-up table that stores for each number of iterations the
  // corresponding inlier ratio.
  double log_failure_prob = std::log(
      parameters.ransac_parameters_.failure_probability);
  uint64_t max_iterations = GetMaxRANSACIterations(
      parameters.ransac_parameters_.min_inlier_ratio, log_failure_prob,
      parameters.ransac_parameters_.use_T_1_1_test);
  inlier_ratio_per_num_iterations_.resize(max_iterations + 1);
  inlier_ratio_per_num_iterations_[0] = 1.0;
  for (uint64_t i = 1; i <= max_iterations; ++i) {
    double failure_prob_root = std::pow(
        parameters.ransac_parameters_.failure_probability,
        1.0 / static_cast<double>(i));
    const double pow_val =
        (parameters.ransac_parameters_.use_T_1_1_test ? 4.0 : 3.0);
    inlier_ratio_per_num_iterations_[i] =
        std::pow(1.0 - failure_prob_root, 1.0 / pow_val);
  }
}

void P3Pf::ComputePose(
    const std::vector<Vector2d> &points2D,
    const std::vector<Vector3d> &points3D,
    PnPResult *result) const 
{
  result->clear();
  int num_correspondences = static_cast<int>(points2D.size());

  if (num_correspondences <= 0) return;

  std::mt19937 rand_num_gen;
  if (params_.ransac_parameters_.random_seed >= 0) {
    rand_num_gen.seed(params_.ransac_parameters_.random_seed);
  }

  std::uniform_int_distribution<int> uniform_distribution_matches(
      0, num_correspondences - 1);

  std::uniform_real_distribution<double> uniform_distribution_focal(0.0, 1.0);

  int num_focal_values = static_cast<int>(params_.focal_length_values_.size());

  // The number of samples generated for each focal length.
  std::vector<uint64_t> num_samples_per_focal_length(num_focal_values, 0);
  // The pose solver to use.
  P3PSolver p3p_solver;
  const double squared_thres = params_.ransac_parameters_.squared_inlier_threshold;

  // Initializes the probabilities according to the prior distribution.
  std::vector<double> probabilities;
  probabilities.assign(params_.prior_probabilities_.begin(),
                       params_.prior_probabilities_.end());

  CameraPose temp_pose;

  std::vector<Vector2d> sample_points_2D(3);
  std::vector<Vector3d> sample_points_3D(3);

  uint64_t t = 0;
  result->num_inliers_ = 0;
  double epsilon_best = params_.ransac_parameters_.min_inlier_ratio;

  // The index of the focal length value for which the best model found so far
  // was computed.
  int best_index = -1;
  //save the probabilities update process
  /*std::ofstream os( "runtime_record.txt", std::ios::trunc );
  if ( !os.is_open() ){ 
    LOG( INFO ) << " open run-time record file fail: " << std::endl;
    return;
  }*/
  while (true) 
  {
    ++t;

    // Randomly select a focal length value.
    double random_val = uniform_distribution_focal(rand_num_gen);
    int sel_focal_length = SampleFocalLength(probabilities, random_val);
    CHECK_NE(sel_focal_length, -1) << "Could not select a valid focal length";

    num_samples_per_focal_length[sel_focal_length] += 1;

    p3p_solver.set_focal_length(params_.focal_length_values_[sel_focal_length]);

    // Randomly select three 2D-3D matches.
    for (int i = 0; i < 3; ++i) {
      int rand_num = uniform_distribution_matches(rand_num_gen);
      sample_points_2D[i] = points2D[rand_num];
      sample_points_3D[i] = points3D[rand_num];
    }

    // Generates poses from the sample.
    p3p_solver.ComputePose(sample_points_2D, sample_points_3D);

    // If the candidate solution passes the T{1,1} test then evaluate all
    // correspondences and update the best pose. PassesT11 test returns true if
    // the T11 test is disabled.
    const int rand_num = uniform_distribution_matches(rand_num_gen);
    if (!params_.ransac_parameters_.use_T_1_1_test ||
      p3p_solver.PassesT11Test( points3D[rand_num], points2D[rand_num], squared_thres ) )
    {
      // Evaluates the poses.
      for (int i = 0; i < num_correspondences; ++i) {
        p3p_solver.EvaluateCorrespondence(
            points3D[i], points2D[i],
            params_.ransac_parameters_.squared_inlier_threshold);
      }

      // Tests whether we have found a new best model.
      int num_inliers = p3p_solver.GetBestCameraPose(&temp_pose);

      // Update best model
      if (num_inliers > result->num_inliers_) {
        result->pose_.InitializePose(temp_pose);
        result->num_inliers_ = num_inliers;
        result->sample_points_2D_ = sample_points_2D;
        result->sample_points_3D_ = sample_points_3D;

        double inlier_ratio = static_cast<double>(result->num_inliers_) /
                              static_cast<double>(num_correspondences);

        epsilon_best = max(epsilon_best, inlier_ratio);
        best_index = sel_focal_length;
      }
    }
    //save the probabilities so that we can visualize the evolvement process
    /*for ( auto prob : probabilities ){ 
      os << prob << " ";
    }
    os << std::endl;*/

    // Update the probabilities.
    bool better_inlier_ratio_possible = true;
    const double max_prob_finding_better_model = UpdateSamplingProbabilities(
          num_samples_per_focal_length,
          best_index,
          epsilon_best,
          &probabilities,
          &better_inlier_ratio_possible);

    /*for ( auto prob : probabilities ){
      os << prob << " ";
      }
      os << std::endl;*/

    if (max_prob_finding_better_model == 0.0) {
      // All probabilities are 0, no need to continue sampling.
      break;
    }
  }

  result->num_generated_random_samples_ = t;
}

/*
bool P3Pf::PassesT11Test(const Vector3d& point3D,
                         const Vector2d& point2D,
                         P3PSolver* p3p_solver) const 
{
  p3p_solver->EvaluateCorrespondence(
      point3D,
      point2D,
      params_.ransac_parameters_.squared_inlier_threshold);

  CameraPose temp_pose;
  int num_inliers = p3p_solver->GetBestCameraPose(&temp_pose);

  // Fail the T{1,1} test if the point was not an inlier.
  if (num_inliers <= 0) {
    return false;
  }

  p3p_solver->SetPose(temp_pose);
  return true;
}*/

double P3Pf::UpdateSamplingProbabilities(
      const std::vector<uint64_t>& num_samples_per_value,
      const int best_focal_length_value,
      const double best_inlier_ratio,
      std::vector<double>* probabilities,
      bool* better_inlier_ratio_possible) const {
  double max_prob_finding_better_model;
  if (best_focal_length_value == -1 ||
      best_inlier_ratio <= params_.ransac_parameters_.min_inlier_ratio) {
    max_prob_finding_better_model = ComputeSamplingProbabilitiesIndependent(
        num_samples_per_value,
        best_inlier_ratio,
        probabilities,
        better_inlier_ratio_possible);
  } else {
    max_prob_finding_better_model = ComputeSamplingProbabilitiesDependent(
        num_samples_per_value,
        best_focal_length_value,
        best_inlier_ratio,
        probabilities,
        better_inlier_ratio_possible);
    }
  return max_prob_finding_better_model;
}

double P3Pf::ComputeSamplingProbabilitiesIndependent(
    const std::vector<uint64_t> &num_samples_per_value,
    const double best_inlier_ratio, std::vector<double> *probabilities,
    bool *better_inlier_ratio_possible) const 
{
  int num_focal_length_values = static_cast<int>(params_.focal_length_values_.size());

  double sum_probabilities = 0.0;
  double cdf_best_inliers = QuantizeValue(
      params_.cdf_inlier_ratios_, best_inlier_ratio,
      scaling_factor_inlier_ratios_);

  double max_prob_better_model = 0.0;
  *better_inlier_ratio_possible = false;

  for (int i = 0; i < num_focal_length_values; ++i) {
    // Computes the maximal inlier ratio that can be achieved for the current
    // focal length value.
    double epsilon_max = GetInlierRatioFromNumberIterations(
        num_samples_per_value[i]);
    if (epsilon_max > best_inlier_ratio) {
      // This is required due to using a quantized cdf in order to prevent
      // terminating too early.
      *better_inlier_ratio_possible = true;

      // Computes the probability of finding a better inlier ratio as depending
      // on the cdf of the inlier ratio and multiplies this with the prior
      // probability of the focal length.
      double cdf_epsilon_max = QuantizeValue(
          params_.cdf_inlier_ratios_, epsilon_max,
          scaling_factor_inlier_ratios_);

      double prob_better_model = cdf_epsilon_max - cdf_best_inliers;
      max_prob_better_model = max(max_prob_better_model,
                                       prob_better_model);
      (*probabilities)[i] = prob_better_model * params_.prior_probabilities_[i];
    } else {
      (*probabilities)[i] = 0.0;
    }
    sum_probabilities += (*probabilities)[i];
  }

  if (sum_probabilities > 0.0) {
    for (int i = 0; i < num_focal_length_values; ++i) {
      (*probabilities)[i] /= sum_probabilities;
    }
  }
  return max_prob_better_model;
}

double P3Pf::ComputeSamplingProbabilitiesDependent(
    const std::vector<uint64_t>& num_samples_per_value,
    const int best_focal_length_value,
    const double best_inlier_ratio,
    std::vector<double>* probabilities,
    bool* better_inlier_ratio_possible) const 
{
  int num_focal_length_values = static_cast<int>(params_.focal_length_values_.size());
  CHECK_GE(best_focal_length_value, 0);
  CHECK_LT(best_focal_length_value, num_focal_length_values);

  double sum_probabilities = 0.0;
  double cdf_best_inliers = QuantizeValue(
      params_.cdf_inlier_ratios_, best_inlier_ratio,
      scaling_factor_inlier_ratios_);

  double max_prob_better_model = 0.0;
  *better_inlier_ratio_possible = false;

  // Computes the probability of finding a better model for the currently
  // best focal length value.
  double one_minus_prob_current_model = 1.0;
  uint64_t sum_iterations = 0;
  for (int i = (best_focal_length_value - 1); i >= 0; --i) {
    sum_iterations += num_samples_per_value[i];
    double epsilon_max = GetInlierRatioFromNumberIterations(sum_iterations);

    if (epsilon_max > best_inlier_ratio) {
      *better_inlier_ratio_possible = true;

      // Computes the probability of finding a better inlier ratio as depending
      // on the cdf of the inlier ratio and multiplies this with the prior
      // probability of the focal length.
      double cdf_epsilon_max = QuantizeValue(
          params_.cdf_inlier_ratios_, epsilon_max,
          scaling_factor_inlier_ratios_);

      double prob_better_model = cdf_epsilon_max - cdf_best_inliers;
      max_prob_better_model = max(max_prob_better_model,
                                       prob_better_model);
      (*probabilities)[i] = prob_better_model * params_.prior_probabilities_[i]
          * one_minus_prob_current_model;
    } else {
      (*probabilities)[i] = 0.0;
    }
    sum_probabilities += (*probabilities)[i];
  }

  {
    double epsilon_max = GetInlierRatioFromNumberIterations(
        num_samples_per_value[best_focal_length_value]);
    if (epsilon_max > best_inlier_ratio) {
      // This is required due to using a quantized cdf in order to prevent
      // terminating too early.
      *better_inlier_ratio_possible = true;

      double cdf_epsilon_max = QuantizeValue(
          params_.cdf_inlier_ratios_, epsilon_max,
          scaling_factor_inlier_ratios_);

      double prob_better_model = cdf_epsilon_max - cdf_best_inliers;
      max_prob_better_model = max(max_prob_better_model,
                                       prob_better_model);
      (*probabilities)[best_focal_length_value] = prob_better_model
          * params_.prior_probabilities_[best_focal_length_value];
    } else {
      (*probabilities)[best_focal_length_value] = 0.0;
    }
    sum_probabilities += (*probabilities)[best_focal_length_value];
  }

  sum_iterations = 0;
  for (int i = (best_focal_length_value + 1); i < num_focal_length_values;
      ++i) {
    sum_iterations += num_samples_per_value[i];
    double epsilon_max = GetInlierRatioFromNumberIterations(sum_iterations);
    if (epsilon_max > best_inlier_ratio) {
      *better_inlier_ratio_possible = true;

      // Computes the probability of finding a better inlier ratio as depending
      // on the cdf of the inlier ratio and multiplies this with the prior
      // probability of the focal length.
      double cdf_epsilon_max = QuantizeValue(
          params_.cdf_inlier_ratios_, epsilon_max,
          scaling_factor_inlier_ratios_);

      double prob_better_model = cdf_epsilon_max - cdf_best_inliers;
      max_prob_better_model = max(max_prob_better_model,
                                       prob_better_model);
      (*probabilities)[i] = prob_better_model * params_.prior_probabilities_[i]
          * one_minus_prob_current_model;
    } else {
      (*probabilities)[i] = 0.0;
    }
    sum_probabilities += (*probabilities)[i];
  }

  if (sum_probabilities > 0.0) {
    for (int i = 0; i < num_focal_length_values; ++i) {
      (*probabilities)[i] /= sum_probabilities;
    }
  }

  return max_prob_better_model;
}

int P3Pf::SampleFocalLength(
    const std::vector<double> &probabilities,
    const double random_value) const {
  CHECK_GT(probabilities.size(), 0);
  double sum = 0.0;
  int num_prob = static_cast<int>(probabilities.size());
  for (int i = 0; i < num_prob; ++i) {
    sum += probabilities[i];
    if (random_value < sum) {
      return i;
    }
  }

  CHECK(false) << random_value << " " << sum;

  // This will only be reached if all probabilities are zero, in which case we
  // assume that the focal lengths are independent. Returning -1 will force
  // RANSAC to compute independent probabilities.
  return -1;
}


double P3Pf::GetInlierRatioFromNumberIterations(
    const uint64_t num_iterations) const {
  // This is needed as we may sum up the number of iterations that have been
  // taken for multiple focal length values. In this case, the sum might
  // go over the maximal number of RANSAC iterations defined by
  // params_.ransac_params_.min_inlier_ratio and
  // params_.ransac_params_.failure_probability. Notice that in this case,
  // the corresponding inlier ratio will be at most
  // params_.ransac_params_.min_inlier_ratio. We could return any value here as
  // the probability of finding a better inlier ratio than the minimal inlier
  // ratio is 0.
  if (num_iterations
      > static_cast<uint64_t>(inlier_ratio_per_num_iterations_.size())) {
    return inlier_ratio_per_num_iterations_.back();
  }
  return inlier_ratio_per_num_iterations_[num_iterations];
}

}  // namespace p3pf
