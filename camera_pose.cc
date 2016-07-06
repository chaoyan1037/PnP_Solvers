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

#include "camera_pose.h"

namespace pnp {
using Eigen::Vector3d;
using Eigen::Matrix;
using Eigen::Matrix3d;

CameraPose::CameraPose() {
  R_.setIdentity();
  t_ = Eigen::Vector3d::Zero();
  K_.setIdentity();
  P_.setIdentity();
  T_.setIdentity();
}

CameraPose::~CameraPose() {
}

void CameraPose::InitializePose(const Matrix3d &R,
                                const Vector3d &t,
                                const Matrix3d &K) {
  R_ = R;
  t_ = t;
  K_ = K;

  T_ << R_, t_;
  P_ = K_ * T_;
}

void CameraPose::InitializePose(const CameraPose &pose) {
  InitializePose(pose.R_, pose.t_, pose.K_);
}

}  // namespace p3pf
