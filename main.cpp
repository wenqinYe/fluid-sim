#include <add_force.h>
#include <advect.h>
#include <common.h>
#include <visualization.h>
#include <diffuse.h>
#include <project.h>
#include <dissipate.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <thread>

igl::opengl::glfw::Viewer viz;

// Simulation state
Eigen::VectorXd V_field_x;
Eigen::VectorXd V_field_y;
Eigen::VectorXd V_field_z;

// The scalar field
Eigen::VectorXd S_field;

// User Tuned Parameters
double dt = 0.5; // time step
bool simulating = true;

// Global values also accessible by the functions in src/*
int dim = 30;
bool show_v_field = false;
bool show_s_field = true;

double viscosity = 0.000000001;
double diffusion = 0.000000001;
double dissipation = 0.0;

// Non-User Tuned
double t = 0;         // simulation time
double domain = dim;
int dim3 = std::pow(dim, 3.0);

Eigen::SparseMatrixd laplace_operator(3 * dim3, 3 * dim3);
Eigen::SparseMatrixd laplace_operator_scalar(dim3, dim3);

Eigen::SparseMatrixd divergence_operator(dim3, 3 * dim3);
Eigen::SparseMatrixd gradient_operator(3 * dim3, dim3);

Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> diffusion_solver;
Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> diffusion_solver_scalar;
Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> laplace_solver_scalar;
Eigen::SparseMatrixd diffusion_mat;
Eigen::SparseMatrixd diffusion_mat_scalar;

Eigen::MatrixXd V_global;
Eigen::MatrixXi F_global;
Eigen::MatrixXd C_global(12 * dim3, 4);

int p00;
int p01;
int p02;
int p10;
int p11;
int p12;
int p20;
int p21;
int p22;

void vstep() {
    /******** 1. Apply forces ********/
    Eigen::VectorXd V_field_x_1;
    Eigen::VectorXd V_field_y_1;
    Eigen::VectorXd V_field_z_1;

    Eigen::VectorXd f_x(V_field_x.size());
    Eigen::VectorXd f_y(V_field_y.size());
    Eigen::VectorXd f_z(V_field_z.size());

    f_x.setZero();
    f_y.setZero();
    f_z.setZero();

    for (int z = 1; z < dim/2; z++) {
        f_x(flat_index(dim/2, dim/2, z)) = 0/(t+1.0);
        f_y(flat_index(dim/2, dim/2, z)) = 0/(t+1.0);
        f_z(flat_index(dim/2, dim/2, z)) = 5.0 * std::sin(z) * 1.0/(t+1.0);

        f_x(flat_index(dim/2, dim/2, dim - 1 -z)) = 0/(t+1.0);
        f_y(flat_index(dim/2, dim/2, dim - 1 - z)) = 0/(t+1.0);
        f_z(flat_index(dim/2, dim/2, dim - 1 - z)) = -5.0 * std::sin(z) * 1.0/(t+1.0);
    }

    add_force(V_field_x_1, V_field_x, f_x, dt);
    add_force(V_field_y_1, V_field_y, f_y, dt);
    add_force(V_field_z_1, V_field_z, f_z, dt);

    apply_fixed_boundary_constraint(V_field_x_1, V_field_y_1, V_field_z_1);

    /******** 2. Advect ********/
    Eigen::VectorXd V_field_x_2 = V_field_x_1;
    Eigen::VectorXd V_field_y_2 = V_field_y_1;
    Eigen::VectorXd V_field_z_2 = V_field_z_1;

    advect(
        V_field_x_2, V_field_y_2, V_field_z_2,  // Output vector field
        V_field_x_1, V_field_y_1, V_field_z_1,  // Input vector field
        dt
    );

    /******** 3. Diffuse ********/
    Eigen::VectorXd V_field_x_3;
    Eigen::VectorXd V_field_y_3;
    Eigen::VectorXd V_field_z_3;

    diffuse_scalar(V_field_x_3, V_field_x_2, dt, viscosity);
    diffuse_scalar(V_field_y_3, V_field_y_2, dt, viscosity);
    diffuse_scalar(V_field_z_3, V_field_z_2, dt, viscosity);
    
    apply_fixed_boundary_constraint(V_field_x_3, V_field_y_3, V_field_z_3); 

    /******** 4. Project ********/
    Eigen::VectorXd V_field_x_4;
    Eigen::VectorXd V_field_y_4;
    Eigen::VectorXd V_field_z_4;

    project( 
        V_field_x_4, V_field_y_4, V_field_z_4,  // Output vector field
        V_field_x_3, V_field_y_3, V_field_z_3  // Input vector field
    );

    /******* Update Velocity fields with new values ********/
    V_field_x = V_field_x_4;
    V_field_y = V_field_y_4;
    V_field_z = V_field_z_4;
}

void sstep() {
    /******** 1. Apply forces ********/
    Eigen::VectorXd S_field_1;

    Eigen::VectorXd f_source(S_field.size());
    f_source.setZero();
    for (int i = dim/2-2; i < dim/2+2; i++) {
        for (int j = dim/2-2; j < dim/2+2; j++) {
            for (int k = 1; k < 3; k++) {
                f_source(flat_index(i, j, k)) = 5.0;
                f_source(flat_index(i, j, dim -1 - k)) = 5.0;
            }
        }
    }

    add_force(S_field_1, S_field, f_source, dt);

    apply_fixed_boundary_constraint_scalar(S_field_1);

    /******** 2. Advect ********/
    Eigen::VectorXd S_field_2;

    advect_scalar(
        S_field_2, // Output vector field
        S_field_1,  // Input vector field
        V_field_x, V_field_y, V_field_z,
        dt
    );

    /******** 3. Diffuse ********/
    Eigen::VectorXd S_field_3;

    diffuse_scalar(S_field_3, S_field_2, dt, diffusion);
    
    apply_fixed_boundary_constraint_scalar(S_field_3); 

    /******** 4. Dissipate ********/
    Eigen::VectorXd S_field_4;
    dissipate(S_field_4, S_field_3, dt, dissipation);

    apply_fixed_boundary_constraint_scalar(S_field_4); 

    /******* Update Scalar field with new values ********/
    S_field = S_field_4;
}

void draw_vector_field() {
    for (int k = 0; k < dim; k++) {  // Put k on outside to optimize memory access of V_field
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                int flat = flat_index(i, j, k);

                Eigen::RowVector3d V_magnitude_vector(V_field_x(flat), V_field_y(flat), V_field_z(flat));
                V_magnitude_vector = V_magnitude_vector / V_magnitude_vector.norm();
                V_magnitude_vector *= 0.25;

                Eigen::RowVector3d V1((domain / (double)dim) * (double)i, (domain / (double)dim) * (double)j, (domain / (double)dim) * (double)k);
                Eigen::RowVector3d V2 = V1 + V_magnitude_vector;

                int half_dim = dim / 2;
                Eigen::RowVector3d offset;
                offset << half_dim, half_dim, half_dim;

                if (flat == p00 || flat == p01 || flat == p02 || flat == p10 || flat == p11 || flat == p12 || flat == p20 || flat == p21 || flat == p22) {
                    viz.data().add_edges(V1 - offset, V2 - offset, Eigen::RowVector3d(255, 0, 0));
                    viz.data().point_size = 5;
                }
                else {
                    viz.data().add_edges(V1 - offset, V2 - offset, Eigen::RowVector3d(0, 0, 0));
                    viz.data().point_size = 5;
                }
            }
        }
    }
}


void draw_scalar_field(igl::opengl::glfw::Viewer &viewer) {
    Eigen::MatrixXd C_global(12 * dim3, 4);

    int curr = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                Eigen::RowVector3d center((domain / (double)dim) * (double)i, (domain / (double)dim) * (double)j, (domain / (double)dim) * (double)k);
                double voxel_dim = domain/(double)dim;

                double scalar_val = S_field(flat_index(i, j, k));
                Eigen::RowVectorXd color(4);
                color << 1. , 1. , 1. , scalar_val * 0.25;

                C_global.block<12, 4>(12 * curr, 0) = color.replicate(12, 1);
                curr += 1;
            }
        }
    }
    
    viewer.core().lighting_factor = 0;
    viewer.data().set_colors(C_global);
}

bool simulation_callback() {

    while (simulating) {
        vstep();
        if (show_s_field)
            sstep();
        t += dt;
    }

    return false;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer) {
    if (show_v_field)
        draw_vector_field();
    if (show_s_field)
        draw_scalar_field(viewer);


    return false;
}

void draw_initial_scalar_field() {
   Eigen::MatrixXd black(dim3, 4);
    Eigen::RowVectorXd color(4);
    color << 0. , 0. , 0. , 0.; 


    black = color.replicate(dim3, 1);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/cube.obj", V, F);

    V_global = V.replicate(dim3, 1);
    F_global = F.replicate(dim3, 1);

    int rows = V.rows();
    int cols = V.cols();

    int curr = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                Eigen::RowVector3d center((domain / (double)dim) * (double)i, (domain / (double)dim) * (double)j, (domain / (double)dim) * (double)k);
                double voxel_dim = domain/(double)dim;

                V_global.block<8, 3>(8 * curr, 0) += center.replicate(V.rows(), 1);
                F_global.block<12, 3>(12 * curr, 0).array() += curr * 8;
                C_global.block<12, 4>(12 * curr, 0) = color.replicate(12, 1);
                curr += 1;
            
            }
        }
    }

    viz.data().set_mesh(V_global,F_global);
    viz.data().set_colors(C_global);
    viz.data().show_lines = false;
}

int main(int argc, char **argv) {

    p00 = flat_index((int)(dim / 2)-1, (int)(dim / 2)-1, dim-1);
    p01 = flat_index((int)(dim / 2)-1, (int)(dim / 2), dim-1);
    p02 = flat_index((int)(dim / 2)-1, (int)(dim / 2)+1, dim-1);
    p10 = flat_index((int)(dim / 2)+0, (int)(dim / 2)-1, dim-1);
    p11 = flat_index((int)(dim / 2)+0, (int)(dim / 2), dim-1);
    p12 = flat_index((int)(dim / 2)+0, (int)(dim / 2)+1, dim-1);
    p20 = flat_index((int)(dim / 2)+1, (int)(dim / 2)-1, dim-1);
    p21 = flat_index((int)(dim / 2)+1, (int)(dim / 2), dim-1);
    p22 = flat_index((int)(dim / 2)+1, (int)(dim / 2)+1, dim-1);

    // Construct the laplace operator matrix
    build_laplace_op();
    build_divergence_op();
    build_gradient_op();
    compute_solvers(dt, diffusion);

    viz.core().background_color(0) = 0.529;
    viz.core().background_color(1) = 0.807;
    viz.core().background_color(2) = 0.921;


    /******** Create the draw callback ********/
    viz.callback_post_draw = &draw_callback;

    /******** Add particles to the visualization ********/
    draw_initial_scalar_field();

    /******** Create the initial scalar, velocity field and add them to the visualization ********/
    // Each vector in the velocity field is represented by an entry in each of these vectors (one entry per dimension)
    V_field_x = Eigen::VectorXd(std::pow(dim, 3.0));
    V_field_y = Eigen::VectorXd(std::pow(dim, 3.0));
    V_field_z = Eigen::VectorXd(std::pow(dim, 3.0));

    S_field = Eigen::VectorXd(std::pow(dim, 3.0));
    S_field.setZero();

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                int flat = flat_index(i, j, k);
                V_field_x(flat) = 0;
                V_field_y(flat) = 0;
                V_field_z(flat) = 0;
            }
        }
    }


    apply_fixed_boundary_constraint(
        V_field_x, V_field_y, V_field_z
    );

    // Show the velocity_field in the visualizaiton
    //
    // Indexing of velocity field (i = column, j = row, k = depth/height)
    //   (i -> x axis), (j -> y axis), (k -> z axis):
    //
    // k
    // ^   j
    // |  7
    // | /
    // |/
    // +-----------> i
    if (show_v_field)
        draw_vector_field();


    /******** Start simulation on background thread ********/
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    viz.core().is_animating = true;
    viz.launch();
    return 1;
}
