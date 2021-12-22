#include <EigenTypes.h>
#include <advect.h>
#include <common.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>

/**
 * Applies the advection step on a vector field
 *
 * Input:
 *    v_field_x0 - The x component of the vector field to update
 *    v_field_y0 - The y component of the vector field to update
 *    V_field_z0 - The z component of the vector field to update
 *    dt - The timestep size
 * Output:
 *    v_field_x1 - The x component of the updated vector field
 *    v_field_y1 - The y component of the updated vector field
 *    V_field_z1 - The z component of the updated vector field
 */
void advect(
    Eigen::VectorXd &V_field_x1, Eigen::VectorXd &V_field_y1, Eigen::VectorXd &V_field_z1,  // Output vector field
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0,  // Input vector field
    double dt) {
    V_field_x1 = Eigen::VectorXd(V_field_x0.size());
    V_field_y1 = Eigen::VectorXd(V_field_y0.size());
    V_field_z1 = Eigen::VectorXd(V_field_z0.size());
    double voxel_dim = domain / (double)dim;
    double half_voxel_dim = voxel_dim / 2.0;

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                int index = flat_index(i, j, k);

                // Position represents the center of the voxel
                Eigen::Vector3d position(
                    voxel_dim * (double)i + half_voxel_dim,
                    voxel_dim * (double)j + half_voxel_dim,
                    voxel_dim * (double)k + half_voxel_dim);

                // Backtrace the position some time dt ago using forward euler
                Eigen::Vector3d velocity(V_field_x0(index), V_field_y0(index), V_field_z0(index));
                Eigen::Vector3d backtraced_pos = position - dt * velocity;

                // Interpolate to get new velocity value
                Eigen::Vector3d new_velocity;
                trilinear_interpolation(new_velocity, backtraced_pos, V_field_x0, V_field_y0, V_field_z0);

                // Update velocity using the interpolated value from some time dt ago.
                V_field_x1(index) = new_velocity(0);
                V_field_y1(index) = new_velocity(1);
                V_field_z1(index) = new_velocity(2);
            }
        }
    }

    apply_fixed_boundary_constraint(V_field_x1, V_field_y1, V_field_z1);
}

/**
 * Applies the advection step on a scalar field using a given velocity field.
 *
 * Input:
 *    p_field_in - The scalar field to update
 *    v_field_x0 - The x component of the velocity field
 *    v_field_y0 - The y component of the velocity field
 *    V_field_z0 - The z component of the velocity field
 *    dt - The timestep size
 * Output:
 *    p_field_out - The updated scalar field
 */
void advect_scalar(
    Eigen::VectorXd &p_field_out,  // Output vector field
    Eigen::VectorXd &p_field_in,
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0,  // Input vector field
    double dt) {
    p_field_out = Eigen::VectorXd(p_field_in.size());
    double voxel_dim = domain / (double)dim;
    double half_voxel_dim = voxel_dim / 2.0;

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                int index = flat_index(i, j, k);

                // Position represents the center of the voxel
                Eigen::Vector3d position(
                    voxel_dim * (double)i + half_voxel_dim,
                    voxel_dim * (double)j + half_voxel_dim,
                    voxel_dim * (double)k + half_voxel_dim);

                // Backtrace position using forwrd euler
                Eigen::Vector3d velocity(V_field_x0(index), V_field_y0(index), V_field_z0(index));
                Eigen::Vector3d backtraced_pos = position - dt * velocity;

                // Interpolate to get new velocity value
                double new_scalar;
                trilinear_interpolation_scalar(new_scalar, backtraced_pos, p_field_in);

                p_field_out(index) = new_scalar;
            }
        }
    }
    apply_fixed_boundary_constraint_scalar(p_field_out);
}