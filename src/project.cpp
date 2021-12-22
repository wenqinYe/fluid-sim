#include <EigenTypes.h>
#include <common.h>
#include <project.h>

#include <Eigen/Dense>
#include <iostream>

/**
 * Applies the projection step on a vector field
 *
 * Input:
 *    v_field_x0 - The x component of the vector field to update
 *    v_field_y0 - The y component of the vector field to update
 *    V_field_z0 - The z component of the vector field to update
 * Output:
 *    v_field_x1 - The x component of the updated vector field
 *    v_field_y1 - The y component of the updated vector field
 *    V_field_z1 - The z component of the updated vector field
 */
void project(
    Eigen::VectorXd &V_field_x1, Eigen::VectorXd &V_field_y1, Eigen::VectorXd &V_field_z1,  // Output vector field
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0   // Input vector field
) {
    // This function assumes that we concat V_field_x0, V_field_y0, V_field_z0,
    // into a single 3 * dim3 vector (recall each V_field is of size dim3)
    Eigen::VectorXd global_field(3 * dim3);
    global_field << V_field_x0, V_field_y0, V_field_z0;

    Eigen::VectorXd divW3 = divergence_operator * global_field;

    Eigen::VectorXd q(dim3);
    q = laplace_solver_scalar.solve(divW3);

    if (laplace_solver_scalar.info() != Eigen::Success) {
        std::cout << "solving failed in projection" << std::endl;
        return;
    }

    Eigen::VectorXd new_global_field(3 * dim3);
    new_global_field = global_field - (gradient_operator * q);

    // Update the velocity field with the new velocity fields
    V_field_x1 = new_global_field.segment(0, dim3);
    V_field_y1 = new_global_field.segment(dim3, dim3);
    V_field_z1 = new_global_field.segment(2 * dim3, dim3);

    apply_fixed_boundary_constraint(V_field_x1, V_field_y1, V_field_z1);
}