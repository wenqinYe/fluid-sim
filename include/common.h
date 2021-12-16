#include <Eigen/Dense>
#include <EigenTypes.h>

extern int dim;
extern int dim3;
extern double domain;

int flat_index(int i, int j, int k); 

void apply_fixed_boundary_constraint(
    Eigen::VectorXd &V_field_x, 
    Eigen::VectorXd &V_field_y, 
    Eigen::VectorXd &V_field_z
); 

void trilinear_interpolation(
    Eigen::Vector3d &result,
    Eigen::Vector3d &position, // not index positions, but the actual position in 3d space 
    Eigen::VectorXd &V_field_x0, 
    Eigen::VectorXd &V_field_y0, 
    Eigen::VectorXd &V_field_z0
);
