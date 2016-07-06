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

#ifndef P3PF_SRC_CAMERA_POSE_H_
#define P3PF_SRC_CAMERA_POSE_H_

#include <Eigen/Dense>
#include <vector>
#include "eigen_helpers.h"

namespace pnp {
using Eigen::Vector4d;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;
using Eigen::Matrix;


// Class modeling the pose of a camera, including the internal calibration of
// the camera. The internal calibration is described by the 3x3 matrix K and
// up to four radial distortion parameters. The type of distortion model, i.e.,
// the normal polynomial one or the division distortion model, is not defined
// for the pose. Instead, each application using this class is responsible for
// applying the parameters in a meaningful way.
class CameraPose {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  CameraPose();
  ~CameraPose();

  // Sets the pose of the camera to be K[R|t], following the notation commonly
  // used in Computer Vision.
  void InitializePose(const Matrix3d& R, const Vector3d& t, const Matrix3d& K);

  void InitializePose(const CameraPose &pose);

  void rotation_matrix(Matrix3d *R) const {
    *R = R_;
  }

  void translation(Vector3d *t) const {
    *t = t_;
  }

  void internal_calibration(Matrix3d *K) const {
    *K = K_;
  }

  double focal_length() const {
    return K_(0, 0) / K_(2, 2);
  }


  // Returns the projection matrix P=K[R|t] defined by the pose.
  void projection_matrix(Eigen::Matrix<double, 3, 4> *P) const {
    *P = P_;
  }

  // Returns the transformation matrix T=[R|t], transforming points from the
  // world coordinate system into the local camera coordinate system.
  void transformation_matrix(Eigen::Matrix<double, 3, 4> *T) const {
    *T = T_;
  }

  // Transforms a 3D point p_c in the local camera coordinate system into a
  // point p_w in the world coordinate system.
  void wTc(const Vector3d &p_c, Vector3d *p_w) const {
    *p_w = R_.transpose() * (p_c - t_);
  }

  // Transforms a 3D point p_w in world coordinates into the point p_c given
  // in camera coordinates.
  void cTw(const Vector3d &p_w, Vector3d *p_c) const {
    *p_c = T_ * p_w.homogeneous();
  }

  // Projects a 3D point p_w in world coordinates into the image using only the
  // projection matrix of the camera and ignoring the radial distortion.
  // Returns true if the point is in front of the camera and false otherwise.
  bool iTw(const Vector3d &p_w, Vector2d *p_i) const {
    const Vector3d p_i_h = P_ * p_w.homogeneous();
    // Return false if behind the camera.
    if ((K_(0, 0) * p_i_h[2]) <= 0.0) {
      return false;
    }
    *p_i = p_i_h.hnormalized();

    return true;
  }

  // Projects a 3D point p_w in world coordinates into the image using only the
  // projection matrix of the camera and ignoring the radial distortion.
  // Returns true if the point is in front of the camera and false otherwise.
  // The second parameter returns the 2D image position in homogeneous
  // coordinates.
  bool iTwHomogeneous(const Vector3d &p_w, Vector3d *p_i_h) const {
    *p_i_h = P_ * p_w.homogeneous();

    // Return false if behind the camera.
    if ((K_(0, 0) * (*p_i_h)[2]) <= 0.0) {
      return false;
    }
    return true;
  }

 protected:
  // The rotation matrix of the pose.
  Eigen::Matrix3d R_;
  // The translational part of the pose.
  Eigen::Vector3d t_;
  // The internal calibration of the camera. K_ has the form diag(f, f, 1),
  // where f is the focal length of the camera.
  Eigen::Matrix3d K_;
  // The projection matrix P=K[R|t].
  Eigen::Matrix<double, 3, 4> P_;
  // The transformation matrix T=[R|t] transforming points from the local
  // into the global camera coordinate system.
  Eigen::Matrix<double, 3, 4> T_;
};

// Structure describing information about a camera, including its extrinsic
// and intrinsic calibration and 2D-3D matches.
struct Camera {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // The pose of the camera.
  CameraPose pose;
  // The width of the corresponding image.
  int width;
  // The heigh of the corresponding image.
  int height;
  // The set of 2D image measurements (feature positions). The image positions
  // might be distorted.
  std::vector<Vector2d> pos_2D;
  // The set of undistorted 2D image measurements (feature positions). The i-th
  // entry in pos_2D_undistorted corresponds to the distorted position
  // pos_2D[i].
  std::vector<Vector2d> pos_2D_undistorted;
  // The 3D points that correspond (matched) to the 2D image measurements.
  // pos_3D[i] is the 3D point that forms a 2D-3D match together with pos_2D[i].
  std::vector<Vector3d> pos_3D;
  // Stores for each matching 3D point its id.
  std::vector<int> pos_3D_ids;
};

}  // namespace p3pf

#endif  // P3PF_SRC_CAMERA_POSE_H_
