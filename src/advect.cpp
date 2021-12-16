#include <advect.h>
#include <Eigen/Dense>
#include <EigenTypes.h>

void trilinear_interpolation() {

}

void advect(
    Eigen::VectorXd &V_field_x1, Eigen::VectorXd &V_field_y1, Eigen::VectorXd &V_field_z1, // Output vector field
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0, // Input vector field
    int dim, double domain
) {
    // for (int k = 0; k < dim; k++) {
    //     for (int i = 0; i < dim; i++) {
    //         for (int j = 0; j < dim; j++) {
    //             double voxel_dim = domain/(double)dim);
                
    //             // Position represents the center of the voxel
    //             Eigen::Vector3d position(
    //                 voxel_dim * (double)i + voxel_dim / 2.0, 
    //                 voxel_dim * (double)j + voxel_dim / 2.0,
    //                 voxel_dim * (double)k + voxel_dim / 2.0);
                
    //             // Use RK2 to get backtraced position
                

    //         }
    //     }
    // }
}