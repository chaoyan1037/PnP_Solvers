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

#include "src/p3p_solver.h"

#include <glog/logging.h>
#include <limits>

#include "eigen_helpers.h"

namespace p3pf {
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;

P3PSolver::P3PSolver() : focal_length_(0.0), num_poses_(0) {
  camera_poses_.clear();
  inlier_counters_.clear();
}

P3PSolver::~P3PSolver() {}

int P3PSolver::ComputePose(
    const std::vector<Vector2d> &points2D,
    const std::vector<Vector3d> &points3D) {
  // CHECK_NE(focal_length_, 0.0);
  if (focal_length_ == 0)
    return 0;

  CHECK_EQ(points2D.size(), 3);
  CHECK_EQ(points3D.size(), 3);

  // Sets the internal camera calibration matrix.
  Matrix3d K;
  K << focal_length_, 0.0, 0.0,
       0.0, focal_length_, 0.0,
       0.0, 0.0, 1.0;

  camera_poses_.clear();
  inlier_counters_.clear();

  // Computes the viewing rays corresponding to the 2D points.
  std::vector<Vector3d> viewing_rays(3);
  for (int i = 0; i < 3; ++i) {
    viewing_rays[i] = Vector3d(points2D[i](0), points2D[i](1), focal_length_);
    viewing_rays[i].normalize();
  }

  // Call P3P solver to compute calibrated camera pose.
  std::vector<Matrix3d> rotation_matrices;
  std::vector<Vector3d> camera_centers;
  p3pf_ext_p3p_kneip::P3PComputePoses(viewing_rays,
                                      points3D,
                                      &rotation_matrices,
                                      &camera_centers);

  CHECK_EQ(rotation_matrices.size(), camera_centers.size());

  num_poses_ = static_cast<int>(rotation_matrices.size());

  if (!num_poses_) return 0;

  camera_poses_.resize(num_poses_);
  inlier_counters_.resize(num_poses_);
  for (int i = 0; i < num_poses_; ++i) {
    inlier_counters_[i] = 0;
    Vector3d t = - rotation_matrices[i] * camera_centers[i];
    camera_poses_[i].InitializePose(rotation_matrices[i], t, K);
  }

  return num_poses_;
}

double P3PSolver::ComputeSquaredReprojectionError(
    const CameraPose &pose,
    const Vector3d &p_w,
    const Vector2d &p_i) const {
  Vector2d p_i_w;
  if (pose.iTw(p_w, &p_i_w)) {
    return (p_i - p_i_w).squaredNorm();
  } else {
    // Point does not project into the image because it is behind the camera.
    return std::numeric_limits<double>::max();
  }
}

void P3PSolver::EvaluateCorrespondence(
    const Vector3d &p_w,
    const Vector2d &p_i,
    const double squared_inlier_threshold) 
{
  for (int i = 0; i < num_poses_; ++i) {
    const double squared_reproj_err =
        ComputeSquaredReprojectionError(camera_poses_[i], p_w, p_i);
    if (squared_reproj_err < squared_inlier_threshold) {
      inlier_counters_[i] += 1;
    }
  }
}

int P3PSolver::GetBestCameraPose(CameraPose *pose) const {
  if (camera_poses_.empty()) return -1;
  int max_index = 0;
  for (int i = 1; i < num_poses_; ++i) {
    if (inlier_counters_[i] > inlier_counters_[max_index]) {
      max_index = i;
    }
  }
  *pose = camera_poses_[max_index];

  return inlier_counters_[max_index];
}

void P3PSolver::SetPose(const CameraPose &pose) {
  inlier_counters_.clear();
  inlier_counters_.push_back(0);
  camera_poses_.clear();
  camera_poses_.push_back(pose);
  num_poses_ = 1;
}

}  // namespace p3pf
