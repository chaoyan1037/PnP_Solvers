#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <algorithm>

#ifndef GLOG_NO_ABBREVIATED_SEVERITIES
#define GLOG_NO_ABBREVIATED_SEVERITIES
#endif
#if _WIN32
#include <windows.h>
#endif
#include <glog/logging.h>

#include "pnp.h"
#include "p6p/p6p.h"

static void Project( const ProjMatrix& proMat, const Vector4d& X, Vector3d& P )
{
  P = proMat*X;
  P = 1.0 / P[3] * P;
}

void GeneratePointsCorrInstance(
  int num_correct_matches, int num_wrong_matches,
  int width, int height,
  double focal_length,
  int seed,
  std::vector<Eigen::Vector2d>& points2D,
  std::vector<Eigen::Vector3d>& points3D )
{
  int num_matches = num_correct_matches + num_wrong_matches;
  points2D.resize( num_matches );
  points3D.resize( num_matches );

  std::vector<Eigen::Vector2d> points2D_temp( num_matches );
  std::vector<Eigen::Vector3d> points3D_temp( num_matches );

  std::default_random_engine rng( seed );
  std::uniform_real_distribution<double> point_distribution_x(
    0.0, static_cast<double>( width - 1.0 ) );
  std::uniform_real_distribution<double> point_distribution_y(
    0.0, static_cast<double>( height - 1.0 ) );
  std::uniform_real_distribution<double> depth_distribution( 1.0, 100.0 );

  // Generates 2D and 3D points by randomly selecting a position in the image
  // and giving the point a random depth.
  for ( int i = 0; i < num_matches; ++i ) {
    points2D_temp[i][0] = point_distribution_x( rng );
    points2D_temp[i][1] = point_distribution_y( rng );

    points3D_temp[i][0] = points2D_temp[i][0] / focal_length;
    points3D_temp[i][1] = points2D_temp[i][1] / focal_length;
    points3D_temp[i][2] = 1.0;
    points3D_temp[i] *= depth_distribution( rng );
  }

  // Generates outlier matches by assigning another random position in the
  // image that is at least 10 pixels away from the correct position.
  for ( int i = num_correct_matches; i < num_matches; ++i ) {
    while ( true ) {
      Eigen::Vector2d new_pos( point_distribution_x( rng ),
        point_distribution_y( rng ) );
      if ( ( new_pos - points2D_temp[i] ).squaredNorm() > 100.0 ) {
        points2D_temp[i] = new_pos;
        break;
      }
    }
  }

  std::normal_distribution<double> noise( 0.0, 1.0 );
  // add noise to the inlier matches
  for ( int i = 0; i < num_correct_matches; ++i ){
    points2D_temp[i][0] += noise( rng );
    points2D_temp[i][0] = std::min( points2D_temp[i][0], static_cast<double>( width - 1.0 ) );
    points2D_temp[i][0] = std::max( points2D_temp[i][0], 0.0 );
    points2D_temp[i][1] += noise( rng );
    points2D_temp[i][1] = std::min( points2D_temp[i][1], static_cast<double>( height - 1.0 ) );
    points2D_temp[i][1] = std::max( points2D_temp[i][1], 0.0 );
  }

  // By now, all outliers are at the end of the points2D_temp and points3D_temp
  // vectors. We now randomly permute them to get a random distribution.
  std::vector<int> indices( num_matches );
  for ( int i = 0; i < num_matches; ++i ) indices[i] = i;

  std::shuffle( indices.begin(), indices.end(), rng );

  for ( int i = 0; i < num_matches; ++i ) {
    points2D[i] = points2D_temp[indices[i]];
    points3D[i] = points3D_temp[indices[i]];
  }

  // Finally, rotate all points by 60 degrees around the y-axis and translate
  // the camera by (3, 6, 12)^T.
  double cos_a = cos( 60.0 * M_PI / 180.0 );
  double sin_a = sin( 60.0 * M_PI / 180.0 );
  Eigen::Matrix3d R, Rz, RT;
  R << cos_a, 0.0, sin_a,
    0.0, 1.0, 0.0,
    -sin_a, 0.0, cos_a;
  Rz << cos_a, -sin_a, 0.0,
    sin_a, cos_a, 0.0,
    0.0, 0.0, 1.0;
  //R = R*Rz;
  Eigen::Vector3d c( 3.0, 6.0, 12.0 );

  Eigen::Vector3d t;
  RT = R.transpose();
  t = -RT*c;

  for ( int i = 0; i < num_matches; ++i ) {
    points3D[i] = R * points3D[i] + c;
  }

  Eigen::Matrix3d K;
  K.setZero();
  K( 0, 0 ) = K( 1, 1 ) = focal_length;
  K( 2, 2 ) = 1.0;

  ProjMatrix projMat;
  for ( int i = 0; i < 3; ++i ){
    for ( int j = 0; j < 3; ++j ){
      projMat( i, j ) = RT( i, j );
    }
    projMat( i, 3 ) = t( i );
  }
  std::cout << "extrinsic [R t]: " << std::endl << projMat << std::endl;
  projMat = K*projMat;
  std::cout << "projection matrix: " << std::endl << projMat << std::endl;
}

using pnp::CameraPose;
using pnp::PnPMethod;

int main( int argc, char *argv[] )
{
  google::ParseCommandLineFlags( &argc, &argv, true );
  FLAGS_log_dir = "./log"; //there should be a log Folder 
  google::InitGoogleLogging( argv[0] );

  // Initializes the parameters of the pnp algorithm.
  pnp::PnPParameters pnp_params;
  pnp_params.InitPnPMethod( PnPMethod::P6P );

  auto & p6p_params = *( pnp_params.ptr_pnp_params );

  p6p_params.ransac_parameters_.failure_probability = 0.01;
  p6p_params.ransac_parameters_.min_inlier_ratio = 0.1;
  p6p_params.ransac_parameters_.squared_inlier_threshold = 10.0;
  p6p_params.ransac_parameters_.random_seed = 0;
  p6p_params.ransac_parameters_.use_T_1_1_test = true;

  // Creates a random problem instance with an inlier ratio of 0.3 for an
  // image of dimensions 640x480 and a focal length corresponding to an
  // opening angle of 55.5 degrees.
  double focal_length = 320.0 / tan( 55.5 / 2.0 * M_PI / 180.0 );
  std::vector<Eigen::Vector2d> points2D;
  std::vector<Eigen::Vector3d> points3D;
  GeneratePointsCorrInstance( 30, 70, 640, 320, focal_length, 0, points2D, points3D );

  pnp_params.ptr_pnp_params->refine_pose = true;

  pnp::PnPSolver pnp;
  pnp.Init( pnp_params );
  LOG( INFO ) << "Running p6p ";
  pnp::PnPResult result;
  std::chrono::high_resolution_clock::time_point t1 =
    std::chrono::high_resolution_clock::now();
  pnp.ComputePose( points2D, points3D, result );
  std::chrono::high_resolution_clock::time_point t2 =
    std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> run_time = std::chrono::duration_cast<
    std::chrono::duration<double >> ( t2 - t1 );

  LOG( INFO ) << "Expected number of inliers: 30, found " << result.num_inliers_;

  LOG( INFO ) << "Expected focal length: " << focal_length << ", estimated "
    << "focal length: " << result.pose_.focal_length();

  Eigen::Vector3d t;
  Eigen::Matrix3d R;
  result.pose_.rotation_matrix( &R );
  result.pose_.translation( &t );
  Eigen::Vector3d c;
  c = -R.transpose() * t;
  LOG( INFO ) << "Expected camera center: [3, 6, 12], computed camera center: "
    << c.transpose();
  LOG( INFO ) << "Camera center error: " << ( Eigen::Vector3d( 3.0, 6.0, 12.0 ) -
    c ).norm();
  LOG( INFO ) << "P3P(f) took " << result.num_generated_random_samples_
    << " iterations to compute the pose";
  LOG( INFO ) << "Running P3P(f) took " << run_time.count() << " seconds ";
  return 1;
}