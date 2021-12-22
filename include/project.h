#include <EigenTypes.h>

#include <Eigen/Dense>

void project(
    Eigen::VectorXd &V_field_x1, Eigen::VectorXd &V_field_y1, Eigen::VectorXd &V_field_z1,  // Output vector field
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0   // Input vector field
);