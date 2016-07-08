
#ifndef _P6P_H_
#define _P6P_H_

#include <vector>
#include <Eigen/Dense>

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;

using ProjMatrix = Eigen::Matrix<double, 3, 4>;

namespace p6p{

//! Correspondence struct for 2-image correspondences
struct Correspondence 
{
  //! Initializing constructor
  Correspondence( const Vector2d & _p2d, const Vector3d & _p3d )
    : m_point2D( _p2d ), m_point3D( _p3d )
  {}

  //! Feature in first image
  Vector2d m_point2D;

  //! Feature in second image
  Vector3d m_point3D;
};

//! Alias for a vector of correspondences
using CorrespondenceVec = std::vector< Correspondence >;

class P6PSolver{

public:

  //! Default constructor
  P6PSolver();

  //! Destructor
  ~P6PSolver();

  //! Clear all point correspondences
  void clear( void );

  //! Add a point correspondence
  void addCorrespondence( const Vector2d & _point2D,
    const Vector3d & _point3D );

  //! Add a correspondence
  void addCorrespondence( const Correspondence & _corres );

  //! Compute projection matrix using SVD
  bool computeLinear( void );

  //! Compute projection matrix with extended linear matrix
  bool computeLinearNew( void );

  //! Get computed projection matrix
  void getProjectionMatrix( ProjMatrix & _mat )
  { _mat = m_projectionMatrix; }

  //! Set the projection matrix
  void setProjectionMatrix( ProjMatrix & _mat )
  { m_projectionMatrix = _mat; }

  //! Evaluate a point correspondence (compute squared reprojection error)
  double evaluateCorrespondence( const Vector2d & _point2D,
    const Vector3d & _point3D );

  //! get the position and orientation of the camera as two 3D vectors
  bool getPositionAndOrientation( Vector3d &position, Vector3d &orientation );

  //! get the camera position in the WCF.
  void getCameraPosition( Vector3d & position ) const;

  //! decompose the left 3X3 part(K*R) of projection matrix
  //! based on QR decomposition
  bool decomposeProjMatrix( Matrix3d& K, Matrix3d& R) const;

  //! based on Givens rotation matrix
  bool decomposeProjMatrixGivens( Matrix3d& K, Matrix3d& R) const;
  
  //! Adjust the K,R so that K's diagonal element is positive
  //! and det(R) is positive 
  bool AdjustDecomposedKR( Matrix3d& K, Matrix3d&R )const;

protected:

  //! Compute projection matrix using SVD without rescaling the resulting projection matrix
  bool computeLinearNoRescaling( void );

  //! Function that scales original features, returns false if one of the scaling factors is too close to 0
  bool scaleCorrespondences( void );


  //! Vector of all point correspondences
  CorrespondenceVec mv_corresp;

  //! Vector of all point correspondences, scaled such that x \in [-1.0, 1.0]
  CorrespondenceVec mv_correspScaled;

  //! Computed projection matrix
  ProjMatrix m_projectionMatrix;

  //! Scaling matrix for image space
  Matrix3d m_matScaleImgInv;

  //! Scaling matrix for world space
  Matrix4d m_matScaleWorld;

};

}


#endif // !_P6P_H_
