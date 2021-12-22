#include <dissipate.h>

/**
 * Dissipates a scalar field by one timestep
 * 
 * Input:
 *    sclar_field_in - The scalar field to apply dissipation to
 *    dt - The timestep size
 *    dissipation - The dissipation rate
 * Output:
 *    scalar_field_out - The updated scalar field
 */
void dissipate(Eigen::VectorXd &scalar_field_out, Eigen::VectorXd &scalar_field_in, double dt, double dissipation) {
    scalar_field_out = scalar_field_in / (1.0 + dt * dissipation);
}