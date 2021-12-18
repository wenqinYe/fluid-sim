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

    // Construct the gradient operator matrix
    // ! TODO: This should probably be called the divergence and the other one the laplacian
    // ! https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Vector_Calculus_(Corral)/04%3A_Line_and_Surface_Integrals/4.06%3A_Gradient_Divergence_Curl_and_Laplacian
    Eigen::SparseMatrixd grad_operator(3 * dim3, 3 * dim3);
    typedef Eigen::Triplet<double> TRI;
    std::vector<TRI> tripletList_grad;


    double pos_grad_coeff = 1.0 / (2.0 * dim);
    double neg_grad_coeff = -1.0 / (2.0 * dim);

    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                // Add x gradient operator for i, j, k
                int row_ind = flat_index(i, j, k);
                tripletList_grad.push_back(TRI(row_ind, flat_index(i-1, j, k), neg_grad_coeff));
                tripletList_grad.push_back(TRI(row_ind, flat_index(i+1, j, k), pos_grad_coeff));

                // Add y gradient opeartor for i, j, k
                tripletList_grad.push_back(TRI(row_ind + dim, flat_index(i, j-1, k), neg_grad_coeff)); // ! TODO: the second parameter should be shifted
                tripletList_grad.push_back(TRI(row_ind + dim, flat_index(i, j+1, k), pos_grad_coeff));

                // Add z gradient operator for i, j, k
                tripletList_grad.push_back(TRI(row_ind + 2 * dim, flat_index(i, j, k-1), neg_grad_coeff));
                tripletList_grad.push_back(TRI(row_ind + 2 * dim, flat_index(i, j, k+1), pos_grad_coeff));
            }
        }
    }
    grad_operator.setFromTriplets(tripletList_grad.begin(), tripletList_grad.end());


    // Construct the divergence operator marix
    // ! TODO: This is probably actually called the Laplacian and is computed differently (See above)
    Eigen::SparseMatrixd div_operator(3 * dim3, 3 * dim3);
    std::vector<TRI> tripletList_div;    
    for (int k = 1; k < dim - 1; k++) {
        for (int i = 1; i < dim - 1; i++) {
            for (int j = 1; j < dim - 1; j++) {
                // Add divergence operator to x, y and z components
                for (int component = 0; component < 3; component++) {
                    int row = flat_index(i, j, k) + component * dim3;

                    // ! TODO: Im not entirly sure what happens here. It looks like its doing a first derivative but that doesnt seem to make sense
                    tripletList_div.push_back(TRI(row, flat_index(i-1, j, k), neg_grad_coeff));
                    tripletList_div.push_back(TRI(row, flat_index(i+1, j, k), pos_grad_coeff));

                    tripletList_div.push_back(TRI(row, flat_index(i, j-1, k), neg_grad_coeff));
                    tripletList_div.push_back(TRI(row, flat_index(i, j+1, k), pos_grad_coeff));

                    tripletList_div.push_back(TRI(row, flat_index(i, j, k-1), neg_grad_coeff));
                    tripletList_div.push_back(TRI(row, flat_index(i, j, k+1), pos_grad_coeff));
                }
            }
        }
    }

    div_operator.setFromTriplets(tripletList_div.begin(), tripletList_div.end());

    // Construct sparse identity matrix
    Eigen::SparseMatrixd identity(3 * dim3, 3 * dim3);
    identity.setIdentity();

    // Create the diffusion operator and then solve for the new velocity field
    Eigen::SparseMatrixd A = identity - viscosity * dt * (div_operator * grad_operator);

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