#include <add_force.h>
#include <cmath>
#include <math.h>
#include <iostream>

/**
 * Applies the add_force step for stable fluid simulation on a scalar field.
 * 
 * To use this function on a vector field simply split the vector field's components
 * into individual scalar fields and then use this function on each scalar field component 
 * 
 * Input:
 *    w0 - DIM3 x 1 Input force field
 *    f - DIM3 x 1 Force field
 * Output:
 *    w1 - DIM3 x 1 Output force field
 */
void add_force(Eigen::VectorXd &w1, Eigen::VectorXd &w0, Eigen::VectorXd f, double dt) {
   w1 = w0 + dt * f;
}