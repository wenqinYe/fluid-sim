#include <EigenTypes.h>
#include <common.h>
#include <diffuse.h>

#include <Eigen/Dense>
#include <iostream>

/**
 * Applies the diffusion step on a vector field.
 *
 * Input:
 *    v_field_x0 - The x component of the vector field to update
 *    v_field_y0 - The y component of the vector field to update
 *    V_field_z0 - The z component of the vector field to update
 *    dt - The timestep size
 *    diffusion - The diffusion parameter
 * Output:
 *    v_field_x1 - The x component of the updated vector field
 *    v_field_y1 - The y component of the updated vector field
 *    V_field_z1 - The z component of the updated vector field
 */
void diffuse(
    Eigen::VectorXd &V_field_x1, Eigen::VectorXd &V_field_y1, Eigen::VectorXd &V_field_z1,  // Output vector field
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0,  // Input vector field
    double dt,
    double diffusion) {
    // This function assumes that we concat V_field_x0, V_field_y0, V_field_z0,
    // into a single 3 * dim3 vector (recall each V_field is of size dim3)
    Eigen::VectorXd global_field(3 * dim3);
    global_field << V_field_x0, V_field_y0, V_field_z0;

    Eigen::VectorXd new_global_field(3 * dim3);
    new_global_field = diffusion_solver.solve(global_field);

    if (diffusion_solver.info() != Eigen::Success) {
        std::cout << "solving failed in diffusion" << std::endl;
        return;
    }

    // Update the velocity field with the new velocity fields
    V_field_x1 = new_global_field.segment(0, dim3);
    V_field_y1 = new_global_field.segment(dim3, dim3);
    V_field_z1 = new_global_field.segment(2 * dim3, dim3);

    apply_fixed_boundary_constraint(V_field_x1, V_field_y1, V_field_z1);
}

/**
 * Applies the diffusion step on a scalar field
 *
 * Input:
 *    w0 - DIM3 x 1 Input force field
 *    diffusion - The diffusion parameter
 *    dt - The timestep size
 * Output:
 *    w1 - DIM3 x 1 Output force field
 */
void diffuse_scalar(Eigen::VectorXd &w1, Eigen::VectorXd &w0, double dt, double diffusion) {
    w1 = diffusion_solver_scalar.solve(w0);

    if (diffusion_solver.info() != Eigen::Success) {
        std::cout << "solving failed in diffusion" << std::endl;
        return;
    }
}