#include <EigenTypes.h>

#include <Eigen/Dense>

extern int dim;
extern int dim3;
extern double domain;
extern double viscosity;

extern Eigen::SparseMatrixd laplace_operator;
extern Eigen::SparseMatrixd laplace_operator_scalar;
extern Eigen::SparseMatrixd divergence_operator;
extern Eigen::SparseMatrixd gradient_operator;

extern Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> diffusion_solver;
extern Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> diffusion_solver_scalar;
extern Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> laplace_solver_scalar;
extern Eigen::SparseMatrixd diffusion_mat;
extern Eigen::SparseMatrixd diffusion_mat_scalar;

int flat_index(int i, int j, int k);

void apply_fixed_boundary_constraint(
    Eigen::VectorXd &V_field_x,
    Eigen::VectorXd &V_field_y,
    Eigen::VectorXd &V_field_z);

void apply_fixed_boundary_constraint_scalar(
    Eigen::VectorXd &P_field);

void trilinear_interpolation(
    Eigen::Vector3d &result,
    Eigen::Vector3d &position,  // not index positions, but the actual position in 3d space
    Eigen::VectorXd &V_field_x0,
    Eigen::VectorXd &V_field_y0,
    Eigen::VectorXd &V_field_z0);

void trilinear_interpolation_scalar(
    double &result,
    Eigen::Vector3d &position,  // not index positions, but the actual position in 3d space
    Eigen::VectorXd &field);

void build_laplace_op();
void build_divergence_op();
void build_gradient_op();
void compute_solvers(double dt, double diffusion);
