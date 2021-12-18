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

    double viscocity = 1.0;

    // Construct the gradient operator matrix
    Eigen::SparseMatrixd grad_operator(3 * dim3, 3 * dim3);
    typedef Eigen::Triplet<double> Temp;
    std::vector<Temp> tripletList_grad;


    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                // Add x gradient operator for i, j, k
                int row = flat_index(i, j, k);
                double delta_x = dim;
                tripletList_grad.push_back(Temp(row, flat_index(i-1, j, k), -1.0 / (2.0 * delta_x)));
                tripletList_grad.push_back(Temp(row, flat_index(i+1, j, k), 1.0 / (2.0 * delta_x)));

                // Add y gradient opeartor for i, j, k
                row = flat_index(i, j, k) + dim3;
                double delta_y = dim;
                tripletList_grad.push_back(Temp(row, flat_index(i, j-1, k), -1.0 / (2.0 * delta_y)));
                tripletList_grad.push_back(Temp(row, flat_index(i, j+1, k), 1.0 / (2.0 * delta_y)));

                // Add z gradient operator for i, j, k
                row = flat_index(i, j, k) + 2 * dim3;
                double delta_z = dim;
                tripletList_grad.push_back(Temp(row, flat_index(i, j, k-1), -1.0 / (2.0 * delta_z)));
                tripletList_grad.push_back(Temp(row, flat_index(i, j, k+1), 1.0 / (2.0 * delta_z)));
            }
        }
    }
    grad_operator.setFromTriplets(tripletList_grad.begin(), tripletList_grad.end());


    // Construct the divergence operator marix
    Eigen::SparseMatrixd div_operator(3 * dim3, 3 * dim3);
    std::vector<Temp> tripletList_div;
    for (int k = 1; k < dim - 1; k++) {

        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                double delta_x = dim;
                double delta_y = dim;
                double delta_z = dim;

                // Add divergence operator to x, y and z components
                for (int component = 0; component < 3; component++) {
                    int row = flat_index(i, j, k) + component * dim3;

                    tripletList_div.push_back(Temp(row, flat_index(i-1, j, k), -1.0 / (2.0 * delta_x)));
                    tripletList_div.push_back(Temp(row, flat_index(i+1, j, k), 1.0 / (2.0 * delta_x)));

                    tripletList_div.push_back(Temp(row, flat_index(i, j-1, k), -1.0 / (2.0 * delta_y)));
                    tripletList_div.push_back(Temp(row, flat_index(i, j+1, k), 1.0 / (2.0 * delta_y)));

                    tripletList_div.push_back(Temp(row, flat_index(i, j, k-1), -1.0 / (2.0 * delta_z)));
                    tripletList_div.push_back(Temp(row, flat_index(i, j, k+1), 1.0 / (2.0 * delta_z)));
                }

            }
        }
    }

    div_operator.setFromTriplets(tripletList_div.begin(), tripletList_div.end());

    // Construct sparse identity matrix
    Eigen::SparseMatrixd identity(3 * dim3, 3 * dim3);
    std::vector<Temp> tripletList_identity;
    for (int i = 0; i < dim3 * 3; i++) {
        tripletList_identity.push_back(Temp(i, i, 1.0));
    }
    identity.setFromTriplets(tripletList_identity.begin(), tripletList_identity.end());

    // Create the diffusion operator and then solve for the new velocity field
    Eigen::SparseMatrixd A = identity - viscocity * dt * (div_operator * grad_operator);

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