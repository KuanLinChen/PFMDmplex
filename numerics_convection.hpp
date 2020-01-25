#pragma once
#include "common.hpp" 
#include <petsc.h>
#include "linear_solvers_structure.hpp"
#include "geometry_structure.hpp"
#include "variable_structure.hpp"
using namespace std;

class CCentJST
{
    
private:
  double  Epsilon_2, Epsilon_4 ;
  double  Local_Lambda_i, Local_Lambda_j, MeanLambda ; /*!< \brief Local eingenvalues. */
  double  Param_p, Param_Kappa_2, Param_Kappa_4 ; /*!< \brief Artificial dissipation parameters. */
  double  *Diff_U, *Diff_Lapl ; /*!< \brief Diference of conservative variables and undivided laplacians. */

 //  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
 //  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
 //  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
 //  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
 //  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
 //  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
 //  Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
 //  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
 // , cte_0, cte_1, /*!< \brief Artificial dissipation values. */
 //    ProjGridVel_i, ProjGridVel_j, ProjGridVel;  /*!< \brief Projected grid velocity. */
 //  stretching; /*!< \brief Stretching factor. */
    
public:

  CGeometry *m ;
  CVariable *var ;
  // CSysSolve *s ;
  
  // Vec Gradient[3];

  Vec Und_Lapl ;
  Vec Residue ;
  void Undivided_Laplacian(Vec);

  void ComputeResidual();

  void Init( CGeometry *, CVariable * ) ;
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST();
    
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST(void){};
    
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  //void ComputeResidual();
};
class CUpwTVD
{
    
private:

public:

  CGeometry *m ;
  CVariable *var ;
  string var_name ;
  // CSysSolve *s ;
  // Vec Gradient[3];

  Vec Residue ;
  Vec Gradient[3];
  void Init( CGeometry *, CVariable *, string ) ;
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwTVD();
  void CalculateGraditntLSQ() ;
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwTVD(void){};
    
  //void ComputeResidual();
};