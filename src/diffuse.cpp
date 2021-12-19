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

void diffuse_scalar(Eigen::VectorXd &w1, Eigen::VectorXd &w0, double dt){
    // Construct sparse identity matrix
    Eigen::SparseMatrixd identity(dim3, dim3);
    identity.setIdentity();

    // Create the diffusion operator and then solve for the new velocity field
    Eigen::SparseMatrixd A = identity - viscosity * dt * laplace_operator_scalar;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
    solver.compute(A);

    w1 = solver.solve(w0);

    if(solver.info()!=Eigen::Success) {
        std::cout << "solving failed in diffusion" << std::endl;
        return;
    }
}