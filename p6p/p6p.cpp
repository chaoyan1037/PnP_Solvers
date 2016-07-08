
#include "p6p/p6p.h"

#include <iostream>


using Eigen::MatrixXd;
using Eigen::MatrixBase;

namespace p6p{ 

static void Project( const ProjMatrix& proMat, const Vector3d& X, Vector2d& P )
{
  Eigen::Vector4d X_;
  Eigen::Vector3d P_;
  X_[0] = X[0]; X_[1] = X[1]; X_[2] = X[2]; X_[0] = 1.0;
  
  P_ = proMat*X_;
  P_ = 1.0 / P_[3] * P_;
  P[0] = P_[0]; P[1] = P_[1];
}

P6PSolver::P6PSolver()
{
  clear();
}

P6PSolver::~P6PSolver()
{

}

void P6PSolver::clear( void )
{
  mv_corresp.clear();
  mv_correspScaled.clear();
}


void P6PSolver::addCorrespondence( const Vector2d & _point2D,
  const Vector3d & _point3D )
{
  // add correspondence
  mv_corresp.push_back( Correspondence( _point2D, _point3D ) );
}


void P6PSolver::addCorrespondence( const Correspondence & _corresp )
{
  // add correspondence
  mv_corresp.push_back( _corresp );
}


bool P6PSolver::scaleCorrespondences( void )
{
  CorrespondenceVec::iterator cIter, cEnd;
  cEnd = mv_corresp.end();

  Vector2d pt1;
  Vector3d pt2;

  Vector2d center2D( 0.0, 0.0 );
  Vector3d center3D( 0.0, 0.0, 0.0 );

  //////////////////////////////////////////////
  // determine c.o.g. in images and (projective) world space
  double n = 1.0;

  for ( cIter = mv_corresp.begin(); cIter != cEnd; ++cIter )
  {
    center2D = cIter->m_point2D * ( 1.0 / n ) + center2D * ( ( n - 1.0 ) / n );
    center3D = cIter->m_point3D * ( 1.0 / n ) + center3D * ( ( n - 1.0 ) / n );
    n += 1.0;
  }

  ///////////////////////////////////////////////
  // move features to COG and compute scaling factors

  mv_correspScaled.clear();

  n = 1.0;
  double scale1 = 0.0;
  double scale2 = 0.0;

  for ( cIter = mv_corresp.begin(); cIter != cEnd; ++cIter )
  {
    pt1 = cIter->m_point2D - center2D;
    pt2 = cIter->m_point3D - center3D;

    scale1 = scale1 * ( ( n - 1.0 ) / n ) + pt1.norm() * ( 1.0 / n );
    scale2 = scale2 * ( ( n - 1.0 ) / n ) + pt2.norm() * ( 1.0 / n );

    mv_correspScaled.push_back( Correspondence( pt1, pt2 ) );

    n += 1.0;
  }

  if ( fabs( scale1 ) < 1e-12 || fabs( scale2 ) < 1e-12 )
    return false;

  scale1 = 1.41421 / scale1; // sqrt(2)
  scale2 = 1.73205 / scale2; // sqrt(3)

  //////////////////////////////////////////////
  // scale coords in both images

  cEnd = mv_correspScaled.end();
  for ( cIter = mv_correspScaled.begin(); cIter != cEnd; ++cIter )
  {
    cIter->m_point2D *= scale1;
    cIter->m_point3D *= scale2;
  }


  ////////////////////////////////////////////
  // create scaling matrices that will be multiplied
  // to result matrix

  m_matScaleImgInv.setZero();
  m_matScaleImgInv( 0, 0 ) = 1.0 / scale1;
  m_matScaleImgInv( 1, 1 ) = 1.0 / scale1;
  m_matScaleImgInv( 0, 2 ) = center2D[0];
  m_matScaleImgInv( 1, 2 ) = center2D[1];
  m_matScaleImgInv( 2, 2 ) = 1.0;

  m_matScaleWorld.setZero();
  m_matScaleWorld( 0, 0 ) = scale2;
  m_matScaleWorld( 1, 1 ) = scale2;
  m_matScaleWorld( 2, 2 ) = scale2;
  m_matScaleWorld( 0, 3 ) = -center3D[0] * scale2;
  m_matScaleWorld( 1, 3 ) = -center3D[1] * scale2;
  m_matScaleWorld( 2, 3 ) = -center3D[2] * scale2;
  m_matScaleWorld( 3, 3 ) = 1.0;

  return true;
}


bool P6PSolver::computeLinearNoRescaling( void )
{
  const size_t nrows = 2 * mv_corresp.size();
  const size_t ncols = 12;

  MatrixXd mat_A( nrows, ncols );
  MatrixXd mat_V( ncols, ncols );
  
  mat_A.setZero();
  mat_V.setZero();

  //gmm::clear( mat_A );

  ///////////////////////////////////////////////////
  // coordinate scaling
  if ( !scaleCorrespondences() )
    return false;

  size_t corr = 0, endCorr = mv_correspScaled.size();

  ///////////////////////////////////////////////////
  // generate matrix A, see Hartley & Zisserman, 2nd edition, page 178-179

  for ( corr = 0; corr < endCorr; ++corr )
  {
    mat_A( 2 * corr, 4 ) = -mv_correspScaled[corr].m_point3D[0];
    mat_A( 2 * corr, 5 ) = -mv_correspScaled[corr].m_point3D[1];
    mat_A( 2 * corr, 6 ) = -mv_correspScaled[corr].m_point3D[2];
    mat_A( 2 * corr, 7 ) = -1.0;

    mat_A( 2 * corr, 8 ) = mv_correspScaled[corr].m_point3D[0]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 2 * corr, 9 ) = mv_correspScaled[corr].m_point3D[1]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 2 * corr, 10 ) = mv_correspScaled[corr].m_point3D[2]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 2 * corr, 11 ) = 1.0 * mv_correspScaled[corr].m_point2D[1];


    mat_A( 2 * corr + 1, 0 ) = mv_correspScaled[corr].m_point3D[0];
    mat_A( 2 * corr + 1, 1 ) = mv_correspScaled[corr].m_point3D[1];
    mat_A( 2 * corr + 1, 2 ) = mv_correspScaled[corr].m_point3D[2];
    mat_A( 2 * corr + 1, 3 ) = 1.0;

    mat_A( 2 * corr + 1, 8 ) = -mv_correspScaled[corr].m_point3D[0]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 2 * corr + 1, 9 ) = -mv_correspScaled[corr].m_point3D[1]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 2 * corr + 1, 10 ) = -mv_correspScaled[corr].m_point3D[2]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 2 * corr + 1, 11 ) = -1.0 * mv_correspScaled[corr].m_point2D[0];
  }

  ////////////////////////////////////////////////////
  // solve system A * p = 0
  // simply using the eigenvector belonging
  // to the smallest eigenvalue.
  Eigen::JacobiSVD<MatrixXd> SVD( mat_A, Eigen::ComputeThinV );
  mat_V = SVD.matrixV();

  for ( int i = 0; i < 3; ++i ){
    for ( int j = 0; j < 4; j++ ){
      m_projectionMatrix( i, j ) = mat_V( 4 * i + j, 11 );
    }
  }

  return true;
}



bool P6PSolver::computeLinear( void )
{
  if ( !computeLinearNoRescaling() )
    return false;

  ////////////////////////////////////////////////////////
  // undo coordinate scaling
  m_projectionMatrix = m_matScaleImgInv*m_projectionMatrix*m_matScaleWorld;

  return true;
}



bool P6PSolver::computeLinearNew( void )
{
  const size_t nrows = 3 * mv_corresp.size();
  const size_t ncols = 12;

  MatrixXd mat_A( nrows, ncols );
  MatrixXd mat_V( ncols, ncols );
  mat_A.setZero();
  mat_V.setZero();

  Vector4d vec_point;

  ///////////////////////////////////////////////////
  // coordinate scaling

  if ( !scaleCorrespondences() )
    return false;

  size_t corr, endCorr = mv_correspScaled.size();

  ///////////////////////////////////////////////////
  // generate matrix A

  for ( corr = 0; corr < endCorr; ++corr )
  {
    vec_point[0] = mv_correspScaled[corr].m_point3D[0];
    vec_point[1] = mv_correspScaled[corr].m_point3D[1];
    vec_point[2] = mv_correspScaled[corr].m_point3D[2];
    vec_point[3] = 1.0;

    mat_A( 3 * corr, 0 ) = -vec_point[0];
    mat_A( 3 * corr, 1 ) = -vec_point[1];
    mat_A( 3 * corr, 2 ) = -vec_point[2];
    mat_A( 3 * corr, 3 ) = -vec_point[3];

    mat_A( 3 * corr, 8 ) = vec_point[0]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 3 * corr, 9 ) = vec_point[1]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 3 * corr, 10 ) = vec_point[2]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 3 * corr, 11 ) = vec_point[3]
      * mv_correspScaled[corr].m_point2D[0];


    mat_A( 3 * corr + 1, 4 ) = -vec_point[0];
    mat_A( 3 * corr + 1, 5 ) = -vec_point[1];
    mat_A( 3 * corr + 1, 6 ) = -vec_point[2];
    mat_A( 3 * corr + 1, 7 ) = -vec_point[3];

    mat_A( 3 * corr + 1, 8 ) = vec_point[0]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 3 * corr + 1, 9 ) = vec_point[1]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 3 * corr + 1, 10 ) = vec_point[2]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 3 * corr + 1, 11 ) = vec_point[3]
      * mv_correspScaled[corr].m_point2D[1];



    mat_A( 3 * corr + 2, 0 ) = -vec_point[0]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 3 * corr + 2, 1 ) = -vec_point[1]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 3 * corr + 2, 2 ) = -vec_point[2]
      * mv_correspScaled[corr].m_point2D[1];
    mat_A( 3 * corr + 2, 3 ) = -vec_point[3]
      * mv_correspScaled[corr].m_point2D[1];

    mat_A( 3 * corr + 2, 4 ) = vec_point[0]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 3 * corr + 2, 5 ) = vec_point[1]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 3 * corr + 2, 6 ) = vec_point[2]
      * mv_correspScaled[corr].m_point2D[0];
    mat_A( 3 * corr + 2, 7 ) = vec_point[3]
      * mv_correspScaled[corr].m_point2D[0];
  }

  ////////////////////////////////////////////////////
  // solve system A * p = 0

  // simply using the eigenvector belonging
  // to the smallest eigenvalue.
  Eigen::JacobiSVD<MatrixXd> SVD( mat_A, Eigen::ComputeFullV );
  mat_V = SVD.matrixV();

  // last column of V is solution 
  for ( int i = 0; i < 3; ++i ){
    for ( int j = 0; j < 4; j++ ){
      m_projectionMatrix( i, j ) = mat_V( 4 * i + j, 11 );
    }
  }

  ////////////////////////////////////////////////////////
  // undo coordinate scaling

  m_projectionMatrix = m_matScaleImgInv*m_projectionMatrix*m_matScaleWorld;

  return true;
}



double P6PSolver::evaluateCorrespondence( const Vector2d & _point2D,
  const Vector3d & _point3D )
{
  static double retval, val;
  Vector2d imgPoint;

  Project( m_projectionMatrix, _point3D, imgPoint );

  val = _point2D[0] - imgPoint[0];
  retval = val * val;
  val = _point2D[1] - imgPoint[1];
  retval += val * val;

  return retval;
}


bool P6PSolver::getPositionAndOrientation( Vector3d &position, Vector3d &orientation )
{
  // get the position of the camera as the vector spanning the right null-space of 
  // the projection matrix, see Hartley & Zisserman, 2nd edition, pages 158-159
  Matrix4d mat_V;
  const ProjMatrix& mat_P = m_projectionMatrix;

  Eigen::JacobiSVD<MatrixXd> SVD( mat_P, Eigen::ComputeFullV );
  mat_V = SVD.matrixV();

  // last col of V is camera position
  for ( int i = 0; i < 3; ++i )
    position[i] = mat_V( i, 3 ) / mat_V( 3, 3 );

  // get the viewing direction of the camera, see Hartley & Zisserman, 2nd ed., pages 160 - 161
  // compute the determinant of the 3x3 part of the projection matrix
  double det = mat_P( 0, 0 ) * mat_P( 1, 1 ) * mat_P( 2, 2 ) 
    + mat_P( 0, 1 ) * mat_P( 1, 2 ) * mat_P( 2, 0 ) 
    + mat_P( 0, 2 ) * mat_P( 1, 0 ) * mat_P( 2, 1 ) 
    - mat_P( 0, 2 ) * mat_P( 1, 1 ) * mat_P( 2, 0 )
    - mat_P( 0, 1 ) * mat_P( 1, 0 ) * mat_P( 2, 2 )
    - mat_P( 0, 0 ) * mat_P( 1, 2 ) * mat_P( 2, 1 );

  // remember that the camera in reconstructions computed by Bundler looks 
  // down the negative z-axis instead of the positive z-axis.
  // So we have to multiply the orientation with -1.0
  for ( int i = 0; i < 3; ++i )
    orientation[i] = -mat_P( 2, i ) * det;

  return true;
}

inline Matrix3d& flipud( Matrix3d & mat ){
  std::swap( mat( 0, 0 ), mat( 2, 0 ) );
  std::swap( mat( 0, 1 ), mat( 2, 1 ) );
  std::swap( mat( 0, 2 ), mat( 2, 2 ) );
  return mat;
}

inline Matrix3d& fliplr( Matrix3d & mat ){
  std::swap( mat( 0, 0 ), mat( 0, 2 ) );
  std::swap( mat( 1, 0 ), mat( 1, 2 ) );
  std::swap( mat( 2, 0 ), mat( 2, 2 ) );
  return mat;
}
// the camera coordinate' X and Y axis is same with image
// camera viewing direction is positive Z
bool P6PSolver::decomposeProjMatrix(Matrix3d& K, Matrix3d& R) const
{
  Matrix3d KR = m_projectionMatrix.block<3, 3>( 0, 0 );

  if ( KR.determinant() == 0.0 ){
    std::cerr << " Error, M is singular!" << std::endl;
    return false;
  }

  KR = flipud( KR ).transpose();
  Eigen::HouseholderQR<Matrix3d> QR;
  QR.compute( KR );

  K = QR.matrixQR().triangularView<Eigen::Upper>();
  K.transposeInPlace();
  fliplr( flipud( K ) );

  R = QR.householderQ();
  R.transposeInPlace();
  flipud( R );

  return true;
}

//! get the camera position in the WCF.
void P6PSolver::getCameraPosition( Vector3d & position ) const
{
  Matrix4d mat_V;
  Eigen::JacobiSVD<MatrixXd> SVD( m_projectionMatrix, Eigen::ComputeFullV );
  mat_V = SVD.matrixV();

  for ( int i = 0; i < 3; ++i )
    position[i] = mat_V( i, 3 ) / mat_V( 3, 3 );
}

bool P6PSolver::decomposeProjMatrixGivens( Matrix3d& K, Matrix3d& R ) const
{
  K = m_projectionMatrix.block<3, 3>( 0, 0 );

  if ( K.determinant() == 0.0 ){
    std::cerr << " Error, M is singular!" << std::endl;
    return false;
  }
  Eigen::Matrix3d matRot;
  double s = 0.0, c = 0.0;

  ////////////////////////////////////
  // Cancellation of element 3,2
  c = std::sqrt( K( 2, 1 )*K( 2, 1 ) + K( 2, 2 )*K( 2, 2 ) );
  s = -K( 2, 1 ) / c;
  c = K( 2, 2 ) / c;

  matRot.setZero();
  matRot( 0, 0 ) = 1.0;
  matRot( 1, 1 ) = matRot( 2, 2 ) = c;
  matRot( 2, 1 ) = s;
  matRot( 1, 2 ) = -s;

  K *= matRot;
  R = matRot.transpose();

  ////////////////////////////////////
  // Cancellation of element 3,1
  c = std::sqrt( K( 2, 0 )*K( 2, 0 ) + K( 2, 2 )*K( 2, 2 ) );
  s = K( 2, 0 ) / c;
  c = K( 2, 2 ) / c;

  matRot.setZero();
  matRot( 1, 1 ) = 1.0;
  matRot( 0, 0 ) = matRot( 2, 2 ) = c;
  matRot( 0, 2 ) = s;
  matRot( 2, 0 ) = -s;

  K *= matRot;
  matRot.transposeInPlace();
  R = matRot*R;

  /////////////////////////////////////
  // Cancellation of element 2,1
  c = std::sqrt( K( 1, 0 )*K( 1, 0 ) + K( 1, 1 )*K( 1, 1 ) );
  s = K( 1, 0 ) / c;
  c = -K( 1, 1 ) / c;

  matRot.setZero();
  matRot( 2, 2 ) = 1.0;
  matRot( 0, 0 ) = matRot( 1, 1 ) = c;
  matRot( 1, 0 ) = s;
  matRot( 0, 1 ) = -s;

  K *= matRot;
  matRot.transposeInPlace();
  R = matRot*R;
  
  return true;
}

bool P6PSolver::AdjustDecomposedKR( Matrix3d& K, Matrix3d&R )const
{
  Eigen::Matrix3d T;
  T.setIdentity();
  if ( K( 0, 0 ) < 0.0 ){ T( 0, 0 ) = -1.0; }
  if ( K( 1, 1 ) < 0.0 ){ T( 1, 1 ) = -1.0; }
  if ( K( 2, 2 ) < 0.0 ){ T( 2, 2 ) = -1.0; }
  // make the diagonal of K positive
  K *= T;
  R = T*R;
  // make the determinant of R positive
  if ( R.determinant() < 0.0 ){
    R = -R;
  }

  double scale = 1.0 / K( 2, 2 );
  if ( !std::isinf( scale ) && !std::isnan( scale ) ){
    K *= scale;
    return true;
  }
  else return false;
}

}// namespace p6p

