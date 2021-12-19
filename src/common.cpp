#include <common.h>

int flat_index(int i, int j, int k) {
    return i + (j * dim) + (dim * dim * k); 
}

void apply_fixed_boundary_constraint(
    Eigen::VectorXd &V_field_x, 
    Eigen::VectorXd &V_field_y, 
    Eigen::VectorXd &V_field_z
) {

    // Wall where k is 0
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            V_field_x(flat_index(i, j, 0)) = -1 * V_field_x(flat_index(i, j, 1));
            V_field_y(flat_index(i, j, 0)) = -1 * V_field_y(flat_index(i, j, 1));
            V_field_z(flat_index(i, j, 0)) = -1 * V_field_z(flat_index(i, j, 1));
        }
    }

    // Wall where k is dim - 1
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            V_field_x(flat_index(i, j, dim-1)) = -1 * V_field_x(flat_index(i, j, dim-2));
            V_field_y(flat_index(i, j, dim-1)) = -1 * V_field_y(flat_index(i, j, dim-2));
            V_field_z(flat_index(i, j, dim-1)) = -1 * V_field_z(flat_index(i, j, dim-2));
        }
    }

    // Wall where j is 0
    for (int i = 0; i < dim; i++) {
        for (int k = 0; k < dim; k++) {
            V_field_x(flat_index(i, 0, k)) = -1 * V_field_x(flat_index(i, 1, k));
            V_field_y(flat_index(i, 0, k)) = -1 * V_field_y(flat_index(i, 1, k));
            V_field_z(flat_index(i, 0, k)) = -1 * V_field_z(flat_index(i, 1, k));
        }
    }

    // Wall where j is dim-1
    for (int i = 0; i < dim; i++) {
        for (int k = 0; k < dim; k++) {
            V_field_x(flat_index(i, dim-1, k)) = -1 * V_field_x(flat_index(i, dim-2, k));
            V_field_y(flat_index(i, dim-1, k)) = -1 * V_field_y(flat_index(i, dim-2, k));
            V_field_z(flat_index(i, dim-1, k)) = -1 * V_field_z(flat_index(i, dim-2, k));
        }
    }

    // Wall where i is 0
    for (int j = 0; j < dim; j++) {
        for (int k = 0; k < dim; k++) {
            V_field_x(flat_index(0, j, k)) = -1 * V_field_x(flat_index(1, j, k));
            V_field_y(flat_index(0, j, k)) = -1 * V_field_y(flat_index(1, j, k));
            V_field_z(flat_index(0, j, k)) = -1 * V_field_z(flat_index(1, j, k));
        }
    }

    // Wall where i is dim-1
    for (int j = 0; j < dim; j++) {
        for (int k = 0; k < dim; k++) {
            V_field_x(flat_index(dim-1, j, k)) = -1 * V_field_x(flat_index(dim-2, j, k));
            V_field_y(flat_index(dim-1, j, k)) = -1 * V_field_y(flat_index(dim-2, j, k));
            V_field_z(flat_index(dim-1, j, k)) = -1 * V_field_z(flat_index(dim-2, j, k));
        }
    }

    // The idea here is on the walls, set the field value to be the same as the feild value
    // just outside of that wall (with the relevant normal component zeroed out)
    
    // // Wall where k is 0
    // for (int i = 0; i < dim; i++) {
    //     for (int j = 0; j < dim; j++) {
    //         V_field_x(flat_index(i, j, 0)) = V_field_x(flat_index(i, j, 1));
    //         V_field_y(flat_index(i, j, 0)) = V_field_y(flat_index(i, j, 1));
    //         V_field_z(flat_index(i, j, 0)) = 0;
    //     }
    // }

    // // Wall where k is dim - 1
    // for (int i = 0; i < dim; i++) {
    //     for (int j = 0; j < dim; j++) {
    //         V_field_x(flat_index(i, j, dim-1)) = V_field_x(flat_index(i, j, dim-2));
    //         V_field_y(flat_index(i, j, dim-1)) = V_field_y(flat_index(i, j, dim-2));
    //         V_field_z(flat_index(i, j, dim-1)) = 0;
    //     }
    // }

    // // Wall where j is 0
    // for (int i = 0; i < dim; i++) {
    //     for (int k = 0; k < dim; k++) {
    //         V_field_x(flat_index(i, 0, k)) = V_field_x(flat_index(i, 1, k));
    //         V_field_y(flat_index(i, 0, k)) = 0;
    //         V_field_z(flat_index(i, 0, k)) = V_field_z(flat_index(i, 1, k));
    //     }
    // }

    // // Wall where j is dim-1
    // for (int i = 0; i < dim; i++) {
    //     for (int k = 0; k < dim; k++) {
    //         V_field_x(flat_index(i, dim-1, k)) = V_field_x(flat_index(i, dim-2, k));
    //         V_field_y(flat_index(i, dim-1, k)) = 0;
    //         V_field_z(flat_index(i, dim-1, k)) = V_field_z(flat_index(i, dim-2, k));
    //     }
    // }

    // // Wall where i is 0
    // for (int j = 0; j < dim; j++) {
    //     for (int k = 0; k < dim; k++) {
    //         V_field_x(flat_index(0, j, k)) = 0;
    //         V_field_y(flat_index(0, j, k)) = V_field_y(flat_index(1, j, k));
    //         V_field_z(flat_index(0, j, k)) = V_field_z(flat_index(1, j, k));
    //     }
    // }

    // // Wall where i is dim-1
    // for (int j = 0; j < dim; j++) {
    //     for (int k = 0; k < dim; k++) {
    //         V_field_x(flat_index(dim-1, j, k)) = 0;
    //         V_field_y(flat_index(dim-1, j, k)) = V_field_y(flat_index(dim-2, j, k));
    //         V_field_z(flat_index(dim-1, j, k)) = V_field_z(flat_index(dim-2, j, k));
    //     }
    // }
}

void trilinear_interpolation(
    Eigen::Vector3d &result,
    Eigen::Vector3d &position, // not index positions, but the actual position in 3d space 
    Eigen::VectorXd &V_field_x0, 
    Eigen::VectorXd &V_field_y0, 
    Eigen::VectorXd &V_field_z0
) {
    double voxel_dim = domain/(double)dim;

    // Constrain it so that we are not doing interpolation on the outer edges
    Eigen::Vector3d pos = position;

    // Get the indices for the closest cells
    double x = (pos(0) - voxel_dim/2.0)/voxel_dim;
    double y = (pos(1) - voxel_dim/2.0)/voxel_dim;
    double z = (pos(2) - voxel_dim/2.0)/voxel_dim;

    double max_idx = dim - 2;
    x = std::max(1.0, x);
    x = std::min(max_idx, x);
    y = std::max(1.0, y);
    y = std::min(max_idx, y);
    z = std::max(1.0, z);
    z = std::min(max_idx, z);

    double x_floor = std::floor(x);
    double y_floor = std::floor(y);
    double z_floor = std::floor(z);

    double x_ceil = x_floor + 1;
    double y_ceil = y_floor + 1;
    double z_ceil = z_floor + 1;

    double x_diff = (x - x_floor) / (x_ceil - x_floor);
    double y_diff = (y - y_floor) / (y_ceil - y_floor);
    double z_diff = (z - z_floor) / (z_ceil - z_floor);

    // Get the field values for the 8 cells that are close to the passed in pos
    int idx = flat_index(x_floor, y_floor, z_floor);
    Eigen::Vector3d c000 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    idx = flat_index(x_ceil, y_floor, z_floor);
    Eigen::Vector3d c100 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    idx = flat_index(x_floor, y_ceil, z_floor);
    Eigen::Vector3d c010 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    idx = flat_index(x_ceil, y_ceil, z_floor);
    Eigen::Vector3d c110 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    idx = flat_index(x_floor, y_floor, z_ceil);
    Eigen::Vector3d c001 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    idx = flat_index(x_ceil, y_floor, z_ceil);
    Eigen::Vector3d c101 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    idx = flat_index(x_floor, y_ceil, z_ceil);
    Eigen::Vector3d c011 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    idx = flat_index(x_ceil, y_ceil, z_ceil);
    Eigen::Vector3d c111 = Eigen::Vector3d(V_field_x0(idx), V_field_y0(idx), V_field_z0(idx));

    // Do the trilinear interpolation
    Eigen::Vector3d c00 = c000 * (1-x_diff) + c100 * x_diff;
    Eigen::Vector3d c01 = c001 * (1-x_diff) + c101 * x_diff;
    Eigen::Vector3d c10 = c010 * (1-x_diff) + c110 * x_diff;
    Eigen::Vector3d c11 = c011 * (1-x_diff) + c111 * x_diff;

    Eigen::Vector3d c0 = c00 * (1-y_diff) + c10 * y_diff;
    Eigen::Vector3d c1 = c01 * (1-y_diff) + c11 * y_diff;

    result = c0 * (1-z_diff) + c1 * z_diff; 
}

void build_laplace_op() {
    // Construct the laplace operator matrix
    typedef Eigen::Triplet<double> TRI;
    std::vector<TRI> tripletList_laplace;
    std::vector<TRI> tripletList_laplace_scalar;

    double voxel_dim = domain/(double)dim;

    double pos_grad_coeff = 1.0 / std::pow(2.0 * voxel_dim, 2.0);
    double neg_grad_coeff = -1.0 / std::pow(2.0 * voxel_dim, 2.0);

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                int row_ind = flat_index(i, j, k);
                // Laplacian operator for x component of vector field
                tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j, k), 6 * neg_grad_coeff));
                tripletList_laplace.push_back(TRI(row_ind, flat_index(i-1, j, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(row_ind, flat_index(i+1, j, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j-1, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j+1, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j, k-1), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j, k+1), pos_grad_coeff));

                // Laplacian operator for y component of vector field
                tripletList_laplace.push_back(TRI(dim3 + row_ind, dim3 + flat_index(i, j, k), 6 * neg_grad_coeff));
                tripletList_laplace.push_back(TRI(dim3 + row_ind, dim3 + flat_index(i-1, j, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(dim3 + row_ind, dim3 + flat_index(i+1, j, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(dim3 + row_ind, dim3 + flat_index(i, j-1, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(dim3 + row_ind, dim3 + flat_index(i, j+1, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(dim3 + row_ind, dim3 + flat_index(i, j, k-1), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(dim3 + row_ind, dim3 + flat_index(i, j, k+1), pos_grad_coeff));

                // Laplacian operator for z component of vector field
                tripletList_laplace.push_back(TRI(2 * dim3 + row_ind, 2 * dim3 + flat_index(i, j, k), 6 * neg_grad_coeff));
                tripletList_laplace.push_back(TRI(2 * dim3 + row_ind, 2 * dim3 + flat_index(i-1, j, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(2 * dim3 + row_ind, 2 * dim3 + flat_index(i+1, j, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(2 * dim3 + row_ind, 2 * dim3 + flat_index(i, j-1, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(2 * dim3 + row_ind, 2 * dim3 + flat_index(i, j+1, k), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(2 * dim3 + row_ind, 2 * dim3 + flat_index(i, j, k-1), pos_grad_coeff));
                tripletList_laplace.push_back(TRI(2 * dim3 + row_ind, 2 * dim3 + flat_index(i, j, k+1), pos_grad_coeff));
            
                // Build the scalar version of the laplacian
                tripletList_laplace_scalar.push_back(TRI(row_ind, flat_index(i, j, k), 6 * neg_grad_coeff));
                tripletList_laplace_scalar.push_back(TRI(row_ind, flat_index(i-1, j, k), pos_grad_coeff));
                tripletList_laplace_scalar.push_back(TRI(row_ind, flat_index(i+1, j, k), pos_grad_coeff));
                tripletList_laplace_scalar.push_back(TRI(row_ind, flat_index(i, j-1, k), pos_grad_coeff));
                tripletList_laplace_scalar.push_back(TRI(row_ind, flat_index(i, j+1, k), pos_grad_coeff));
                tripletList_laplace_scalar.push_back(TRI(row_ind, flat_index(i, j, k-1), pos_grad_coeff));
                tripletList_laplace_scalar.push_back(TRI(row_ind, flat_index(i, j, k+1), pos_grad_coeff));
            }
        }
    }

    laplace_operator.setFromTriplets(tripletList_laplace.begin(), tripletList_laplace.end());
    laplace_operator_scalar.setFromTriplets(tripletList_laplace_scalar.begin(), tripletList_laplace_scalar.end());
}

void build_divergence_op() {
    // Construct the divergence operator
    typedef Eigen::Triplet<double> TRI;
    std::vector<TRI> tripletList_divergence;

    double voxel_dim = domain/(double)dim;

    double pos_grad_coeff = 1.0 / (2.0 * voxel_dim);
    double neg_grad_coeff = -1.0 / (2.0 * voxel_dim);

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                int row_ind = flat_index(i, j, k);
                
                tripletList_divergence.push_back(TRI(row_ind, flat_index(i-1, j, k), neg_grad_coeff));
                tripletList_divergence.push_back(TRI(row_ind, flat_index(i+1, j, k), pos_grad_coeff));
                
                tripletList_divergence.push_back(TRI(row_ind, flat_index(i, j-1, k) + dim3, neg_grad_coeff));
                tripletList_divergence.push_back(TRI(row_ind, flat_index(i, j+1, k) + dim3, pos_grad_coeff));

                tripletList_divergence.push_back(TRI(row_ind, flat_index(i, j, k-1) + 2 * dim3, neg_grad_coeff));
                tripletList_divergence.push_back(TRI(row_ind, flat_index(i, j, k+1) + 2 * dim3, pos_grad_coeff));
            }
        }
    }
    divergence_operator.setFromTriplets(tripletList_divergence.begin(), tripletList_divergence.end());
}

void build_gradient_op() {
    // Construct the gradient operator
    typedef Eigen::Triplet<double> TRI;
    std::vector<TRI> tripletList_gradient;

    double voxel_dim = domain/(double)dim;

    double pos_grad_coeff = 1.0 / (2.0 * voxel_dim);
    double neg_grad_coeff = -1.0 / (2.0 * voxel_dim);

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                int row_ind = flat_index(i, j, k);

                tripletList_gradient.push_back(TRI(row_ind, flat_index(i-1, j, k), neg_grad_coeff));
                tripletList_gradient.push_back(TRI(row_ind, flat_index(i+1, j, k), pos_grad_coeff));
                
                tripletList_gradient.push_back(TRI(row_ind + dim3, flat_index(i, j-1, k), neg_grad_coeff));
                tripletList_gradient.push_back(TRI(row_ind + dim3, flat_index(i, j+1, k), pos_grad_coeff));

                tripletList_gradient.push_back(TRI(row_ind + 2 * dim3, flat_index(i, j, k-1), neg_grad_coeff));
                tripletList_gradient.push_back(TRI(row_ind + 2 * dim3, flat_index(i, j, k+1), pos_grad_coeff));
            }
        }
    }
    gradient_operator.setFromTriplets(tripletList_gradient.begin(), tripletList_gradient.end());
}