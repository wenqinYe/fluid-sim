#include <dissipate.h>

void dissipate(Eigen::VectorXd &scalar_field_out, Eigen::VectorXd &scalar_field_in, double dt, double dissipation) {
    scalar_field_out = scalar_field_in / (1.0 + dt * dissipation);
}