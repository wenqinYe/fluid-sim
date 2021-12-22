#include <common.h>

#include <Eigen/Dense>

void dissipate(Eigen::VectorXd &scalar_field_out, Eigen::VectorXd &scalar_field_in, double dt, double dissipation);