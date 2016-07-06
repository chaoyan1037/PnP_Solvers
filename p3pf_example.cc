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

#include <Eigen/Dense>
#include <Eigen/StdVector>

#ifndef GLOG_NO_ABBREVIATED_SEVERITIES
#define GLOG_NO_ABBREVIATED_SEVERITIES
#endif
#if _WIN32
#define NOMINMAX
#include <windows.h>
#endif
#include <glog/logging.h>

#include <algorithm>
#include <chrono>
#include <random>

#include "pnp.h"
#include "pnp_params.h"
#include "camera_pose.h"
#include "eigen_helpers.h"
//#include "p3pf/src/p3pf_params.h"
//#include "p3pf/src/p3pf.h"

using pnp::CameraPose;
using pnp::PnPMethod;

// Creates an instance for P3P(f) by randomly selecting points inside the image,
// associating them with depth values, and transforming the resulting 3D points.
// Also adds a specified number of outliers.
void GenerateP3PfInstance(int num_correct_matches, int num_wrong_matches,
                          int width, int height,
                          double focal_length,
                          int seed,
                          std::vector<Eigen::Vector2d>* points2D,
                          std::vector<Eigen::Vector3d>* points3D) 
{
  int num_matches = num_correct_matches + num_wrong_matches;
  points2D->resize(num_matches);
  points3D->resize(num_matches);

  std::vector<Eigen::Vector2d> points2D_temp(num_matches);
  std::vector<Eigen::Vector3d> points3D_temp(num_matches);

  std::default_random_engine rng(seed);
  std::uniform_real_distribution<double> point_distribution_x(
      0.0, static_cast<double>(width));
  std::uniform_real_distribution<double> point_distribution_y(
      0.0, static_cast<double>(height));
  std::uniform_real_distribution<double> depth_distribution(1.0, 100.0);

  // Generates 2D and 3D points by randomly selecting a position in the image
  // and giving the point a random depth.
  for (int i = 0; i < num_matches; ++i) {
    points2D_temp[i][0] = point_distribution_x(rng);
    points2D_temp[i][1] = point_distribution_y(rng);

    points3D_temp[i][0] = points2D_temp[i][0] / focal_length;
    points3D_temp[i][1] = points2D_temp[i][1] / focal_length;
    points3D_temp[i][2] = 1.0;
    points3D_temp[i] *= depth_distribution(rng);
  }

  // Generates outlier matches by assigning another random position in the
  // image that is at least 10 pixels away from the correct position.
  for (int i = num_correct_matches; i < num_matches; ++i) {
    while (true) {
      Eigen::Vector2d new_pos(point_distribution_x(rng),
                              point_distribution_y(rng));
      if ((new_pos - points2D_temp[i]).squaredNorm() > 100.0) {
        points2D_temp[i] = new_pos;
        break;
      }
    }
  }
  
  // add noise to the inlier matches
  std::normal_distribution<double> noise( 0.0, 1.0 );
  for ( int i = 0; i < num_correct_matches; ++i ){
    points2D_temp[i][0] += noise( rng );
    points2D_temp[i][0] = std::min( points2D_temp[i][0], static_cast<double>( width ) );
    points2D_temp[i][0] = std::max( points2D_temp[i][0], 0.0 );
    points2D_temp[i][1] += noise( rng );
    points2D_temp[i][1] = std::min( points2D_temp[i][1], static_cast<double>( height ) );
    points2D_temp[i][1] = std::max( points2D_temp[i][1], 0.0 );
  }

  // By now, all outliers are at the end of the points2D_temp and points3D_temp
  // vectors. We now randomly permute them to get a random distribution.
  std::vector<int> indices(num_matches);
  for (int i = 0; i < num_matches; ++i) indices[i] = i;

  std::shuffle(indices.begin(), indices.end(), rng);

  for (int i = 0; i < num_matches; ++i) {
    (*points2D)[i] = points2D_temp[indices[i]];
    (*points3D)[i] = points3D_temp[indices[i]];
  }

  // Finally, rotate all points by 60 degrees around the y-axis and translate
  // the camera by (3, 6, 12)^T.
  double cos_a = cos(60.0 * M_PI / 180.0);
  double sin_a = sin(60.0 * M_PI / 180.0);
  Eigen::Matrix3d R;
  R << cos_a, 0.0, sin_a,
       0.0, 1.0, 0.0,
       -sin_a, 0.0, cos_a;
  Eigen::Vector3d c(3.0, 6.0, 12.0);

  for (int i = 0; i < num_matches; ++i) {
    (*points3D)[i] = R * (*points3D)[i] + c;
  }
}

int main(int argc, char *argv[]) 
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_log_dir = "./log"; //there should be a log Folder 
  google::InitGoogleLogging(argv[0]);

  // The original list of opening angles used for the experiments in the
  // paper.
  std::vector<double> opening_angles = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5,
      8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
      20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5,
      32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5,
      44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5, 51.5, 52.5, 53.5, 54.5, 55.5,
      56.5, 57.5, 58.5, 59.5, 60.5, 61.5, 62.5, 63.5, 64.5, 65.5, 66.5, 67.5,
      68.5, 69.5, 70.5, 71.5, 72.5, 73.5, 74.5, 75.5, 76.5, 77.5, 78.5, 79.5,
      80.5, 81.5, 82.5, 83.5, 84.5, 85.5, 86.5, 87.5, 88.5, 89.5, 90.5, 91.5,
      92.5, 93.5, 94.5, 95.5, 96.5, 97.5, 98.5, 99.5};
#if 0 // do not use prior
  std::vector<double> prior_probabilities( 100, 0.01 );
#else
  std::vector<double> prior_probabilities = { 0.0019519, 0.00060173, 0.00044518,
      0.0010176, 0.0048432, 0.0034343, 0.0077197, 0.0038012, 0.0034098,
      0.0062228, 0.0046328, 0.0039968, 0.0050633, 0.0044176, 0.0109, 0.0062325,
      0.0060466, 0.017308, 0.01993, 0.014945, 0.0065994, 0.0080915, 0.015024,
      0.012196, 0.015298, 0.011917, 0.014877, 0.011536, 0.0083802, 0.0090161,
      0.0080915, 0.010332, 0.011599, 0.0089868, 0.00839, 0.0093244, 0.0094515,
      0.0127, 0.0094613, 0.0087569, 0.010621, 0.010763, 0.012118, 0.014549,
      0.013189, 0.014353, 0.01563, 0.02813, 0.050389, 0.053666, 0.047938,
      0.065119, 0.065642, 0.02904, 0.016261, 0.011364, 0.0088645, 0.0071718,
      0.006394, 0.0055574, 0.0067609, 0.010704, 0.022083, 0.025263, 0.015508,
      0.021892, 0.02358, 0.0093928, 0.0039528, 0.0026564, 0.0019813, 0.0021476,
      0.0025586, 0.0023384, 0.001179, 0.0011301, 0.0010322, 0.00092461,
      0.00077295, 0.00072892, 0.00075338, 0.00076806, 0.00068, 0.00089526,
      0.00086101, 0.00082187, 0.00082187, 0.00088058, 0.00091972, 0.00089036,
      0.00079252, 0.00076806, 0.00090015, 0.00088058, 0.0015214, 0.0020987,
      0.0016829, 0.0012328, 0.00091972, 0.00061641 };
#endif
  // Initializes the parameters of the pnp algorithm.
  pnp::PnPParameters pnp_params;
  pnp_params.InitPnPMethod( pnp::PnPMethod::P3PF );

  auto & p3pf_params = *(pnp_params.ptr_p3pf_params);

  p3pf_params.ransac_parameters_.failure_probability = 0.01;
  p3pf_params.ransac_parameters_.min_inlier_ratio = 0.1;
  p3pf_params.ransac_parameters_.squared_inlier_threshold = 10.0;
  p3pf_params.ransac_parameters_.random_seed = 0;
  p3pf_params.ransac_parameters_.use_T_1_1_test = true;

  // Initializes the cdf values. Here, we assume that all inlier ratios have
  // the same probability, i.e., a uniform distribution.
  // Notice that P3P(f) linearly interpolates between the cdf values.
  // If you are using P3P(f) with your own cdf, you need to ensure that the
  // cdf is strictly increasing. Otherwise, P3P(f) might terminate too early.
  // Notice that it is not necessary to specify the inlier ratio values. P3P(f)
  // assumes that the first entry corresponds to an inlier ratio of 0 and the
  // last one to the inlier ratio 1. The intermediate values are assumed to have
  // equal spacing.
  p3pf_params.cdf_inlier_ratios_ = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
    0.8, 0.9, 1.0};

  // Creates a random problem instance with an inlier ratio of 0.3 for an
  // image of dimensions 640x480 and a focal length corresponding to an
  // opening angle of 55.5 degrees.
  double focal_length = 320.0 / tan(55.5 / 2.0 * M_PI / 180.0);
  std::vector<Eigen::Vector2d> points2D;
  std::vector<Eigen::Vector3d> points3D;
  GenerateP3PfInstance(30, 70, 640, 320, focal_length, 0, &points2D, &points3D);

  // Converts the opening angles used for sampling to focal length values.
  p3pf_params.focal_length_values_ = opening_angles;
  for (size_t i = 0; i < opening_angles.size(); ++i) {
    p3pf_params.focal_length_values_[i] = 320.0
        / tan(opening_angles[i] / 2.0 * M_PI / 180.0);
  }

  p3pf_params.prior_probabilities_ = prior_probabilities;

  pnp::PnPSolver pnp;
  pnp.Init( pnp_params );

  // Initializes and runs P3P(f).
  //p3pf::P3Pf p3pf;
  //p3pf.Init(p3pf_params);

  LOG(INFO) << "Running P3P(f) ";
  pnp::PnPResult result;
  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();
  pnp.ComputePose( points2D, points3D, result );
  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> run_time = std::chrono::duration_cast<
      std::chrono::duration<double>>(t2 - t1);

  LOG(INFO) << "Expected number of inliers: 30, found " << result.num_inliers_;

  LOG(INFO) << "Expected focal length: " << focal_length << ", estimated "
            << "focal length: " << result.pose_.focal_length();

  Eigen::Vector3d t;
  Eigen::Matrix3d R;
  result.pose_.rotation_matrix(&R);
  result.pose_.translation(&t);
  Eigen::Vector3d c;
  c = -R.transpose() * t;
  LOG(INFO) << "Expected camera center: [3, 6, 12], computed camera center: "
            << c.transpose();
  LOG(INFO) << "Camera center error: " << (Eigen::Vector3d(3.0, 6.0, 12.0) -
                                             c).norm();
  LOG(INFO) << "P3P(f) took " << result.num_generated_random_samples_
            << " iterations to compute the pose";
  LOG(INFO) << "Running P3P(f) took " << run_time.count() << " seconds ";
}
