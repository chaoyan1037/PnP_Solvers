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

#ifndef P3PF_SRC_P3P_SOLVER_H_
#define P3PF_SRC_P3P_SOLVER_H_

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <vector>

#include "ext/p3p_code_kneip/P3p.h"
#include "camera_pose.h"
#include "eigen_helpers.h"

namespace p3pf {
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;
using pnp::CameraPose;

// Implementation of a PnPSolver class that uses Laurent Kneip's implementation
// of a 3-point pose solver to estimate the external camera parameters from
// three 2D-3D correspondences.
// The 3-point solver is based on the paper
//   L. Kneip, D. Scaramuzza, R. Siegwart. A Novel Parametrization of the
//   Perspective-Three-Point Problem for a Direct Computation of Absolute
//   Camera Position and Orientation. CVPR 2011.
class P3PSolver {
 public:
  P3PSolver();
  ~P3PSolver();

  // Estimates the camera pose from three 2D-3D matches. Returns the number of
  // computed poses or 0 if no pose could be computed. The focal length has to
  // be set using set_focal_length before calling ComputePose. The function
  // returns number of computed poses or 0 if no pose could be computed.
  int ComputePose(
      const std::vector<Vector2d> &points2D,
      const std::vector<Vector3d> &points3D);

  int GetMinimalNumberOfCorrespondences() {
    return 3;
  }

  double focal_length() const {
    return focal_length_;
  }

  // Sets the focal length of the camera. Notice that focal_length can be
  // negative, meaning that the camera is looking down the negative z-axis (as
  // is common in Computer Graphics).
  void set_focal_length(const double focal_length) {
    focal_length_ = focal_length;
  }

  // Evaluates a single correspondence, given by a 3D point p_w in world
  // coordinates and its corresponding 2D image position p_i, against all
  // current poses. If the squared reprojection error is below the given
  // threshold squared_inlier_threshold, the correspondence is counted as an
  // inlier for the pose and the inlier counter of the pose is incremented.
  // This function calls compute_squared_reprojection_error to compute the
  // error.
  void EvaluateCorrespondence(const Vector3d &p_w,
                              const Vector2d &p_i,
                              const double squared_inlier_threshold);

  // Determines and returns the camera pose with the largest number of inliers
  // and this number of inliers. Returns -1 if no camera has been computed.
  int GetBestCameraPose(CameraPose *pose) const;

  // Clears all poses stored in the solver, replaces them with the given pose
  // and resets the inlier counter.
  void SetPose(const CameraPose &pose);

  // Returns the number of poses computed during the last invocation of
  // ComputePose.
  int GetNumberOfCurrentPoses() {
    return num_poses_;
  }

 private:
  // Evaluates a correspondence against the given pose, returns the squared
  // reprojection error. The default implementation uses the iTw function of
  // the camera pose to compute the projection of the 3D point.
  // The parameters specify the camera pose, the position of the 3D point in
  // world coordinates, and the corresponding 2D point position in the image.
  double ComputeSquaredReprojectionError(const CameraPose &pose,
                                         const Vector3d &p_w,
                                         const Vector2d &p_i) const;

  // The focal length used when computing the camera poses.
  double focal_length_;

  // The number of stored poses.
  int num_poses_;

  // Stores the poses computed during the last call of compute_pose.
  std::vector<CameraPose, Eigen::aligned_allocator<CameraPose> > camera_poses_;

  // Stores the inlier counts for the computed poses.
  std::vector<int> inlier_counters_;
};

}  // namespace p3pf

#endif  // P3PF_SRC_P3P_SOLVER_H_
