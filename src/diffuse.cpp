#include <diffuse.h>
#include <Eigen/Dense>
#include <EigenTypes.h>
#include <common.h>
#include <iostream>

void diffuse(
    Eigen::VectorXd &V_field_x1, Eigen::VectorXd &V_field_y1, Eigen::VectorXd &V_field_z1, // Output vector field
    Eigen::VectorXd &V_field_x0, Eigen::VectorXd &V_field_y0, Eigen::VectorXd &V_field_z0, // Input vector field
    double dt
){
    // This function assumes that we concat V_field_x0, V_field_y0, V_field_z0,
    // into a single 3 * dim3 vector (recall each V_field is of size dim3) 
    Eigen::VectorXd global_field(3 * dim3);
    global_field << V_field_x0, V_field_y0, V_field_z0;

    // Construct the laplace pperator matrix
    Eigen::SparseMatrixd laplace_operator(3 * dim3, 3 * dim3);
    typedef Eigen::Triplet<double> TRI;
    std::vector<TRI> tripletList_laplace;

    double pos_grad_coeff = 1.0 / (2.0 * dim);
    double neg_grad_coeff = -1.0 / (2.0 * dim);

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                for (int comp = 0; comp < 3; comp++) {
                    int row_ind = flat_index(i, j, k) + comp * dim3;

                    // Laplacian operator for x component of vector field
                    tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j, k), 6 * neg_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, flat_index(i-1, j, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, flat_index(i+1, j, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j-1, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j+1, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j, k-1), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, flat_index(i, j, k+1), pos_grad_coeff));

                    // Laplacian operator for y component of vector field
                    tripletList_laplace.push_back(TRI(row_ind, dim3 + flat_index(i, j, k), 6 * neg_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, dim3 + flat_index(i-1, j, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, dim3 + flat_index(i+1, j, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, dim3 + flat_index(i, j-1, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, dim3 + flat_index(i, j+1, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, dim3 + flat_index(i, j, k-1), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, dim3 + flat_index(i, j, k+1), pos_grad_coeff));

                    // Laplacian operator for z component of vector field
                    tripletList_laplace.push_back(TRI(row_ind, 2 * dim3 + flat_index(i, j, k), 6 * neg_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, 2 * dim3 + flat_index(i-1, j, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, 2 * dim3 + flat_index(i+1, j, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, 2 * dim3 + flat_index(i, j-1, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, 2 * dim3 + flat_index(i, j+1, k), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, 2 * dim3 + flat_index(i, j, k-1), pos_grad_coeff));
                    tripletList_laplace.push_back(TRI(row_ind, 2 * dim3 + flat_index(i, j, k+1), pos_grad_coeff));
                }
            }
        }
    }
    laplace_operator.setFromTriplets(tripletList_laplace.begin(), tripletList_laplace.end());

    // Construct sparse identity matrix
    Eigen::SparseMatrixd identity(3 * dim3, 3 * dim3);
    identity.setIdentity();

    // Create the diffusion operator and then solve for the new velocity field
    Eigen::SparseMatrixd A = identity - viscosity * dt * laplace_operator;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
    solver.compute(A);

    Eigen::VectorXd new_global_field(3 * dim3);
    new_global_field = solver.solve(global_field);

    if(solver.info()!=Eigen::Success) {
        std::cout << "solving failed in diffusion" << std::endl;
        return;
    }
    
    // Update the velocity field with the new velocity fields
    V_field_x1 = new_global_field.segment(0, dim3);
    V_field_y1 = new_global_field.segment(dim3, dim3);
    V_field_z1 = new_global_field.segment(2 * dim3, dim3);

    apply_fixed_boundary_constraint(V_field_x1, V_field_y1, V_field_z1); 
}