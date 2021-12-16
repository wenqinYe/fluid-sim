#include <advect.h>
#include <Eigen/Dense>
#include <EigenTypes.h>
#include <cmath>
#include <algorithm>
#include <common.h>
#include <iostream>

void advect(
    Eigen::VectorXd &V_field_x1, Eigen::VectorXd &V_field_y1, Eigen::VectorXd &V_field_z1, // Output vector field
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0, // Input vector field
    double dt
) {
    V_field_x1 = Eigen::VectorXd(V_field_x0.size());
    V_field_y1 = Eigen::VectorXd(V_field_y0.size());
    V_field_z1 = Eigen::VectorXd(V_field_z0.size());

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                int index = flat_index(i, j, k);
                double voxel_dim = domain/(double)dim;
                double half_voxel_dim = voxel_dim / 2.0;
                
                // Position represents the center of the voxel
                Eigen::Vector3d position(
                    voxel_dim * (double)i + voxel_dim / 2.0, 
                    voxel_dim * (double)j + voxel_dim / 2.0,
                    voxel_dim * (double)k + voxel_dim / 2.0);
                
                // Start with simple euler (TODO: update to use RK2)
                Eigen::Vector3d velocity(V_field_x0(index), V_field_y0(index), V_field_z0(index));
                Eigen::Vector3d backtraced_pos = position - dt * velocity;

                // Interpolate to get new velocity value
                Eigen::Vector3d new_velocity;
                trilinear_interpolation(new_velocity, backtraced_pos, V_field_x0, V_field_y0, V_field_z0);

                V_field_x1(index) = new_velocity(0);
                V_field_y1(index) = new_velocity(1);
                V_field_z1(index) = new_velocity(2);
            }      
        }
    }

    apply_fixed_boundary_constraint(V_field_x1, V_field_y1, V_field_z1);
}