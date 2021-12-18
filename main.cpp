#include <add_force.h>
#include <advect.h>
#include <common.h>
#include <visualization.h>
#include <diffuse.h>

#include <cmath>
#include <iostream>
#include <thread>

igl::opengl::glfw::Viewer viz;

// Simulation state
Eigen::MatrixXd P;
Eigen::VectorXd V_field_x;
Eigen::VectorXd V_field_y;
Eigen::VectorXd V_field_z;

double t = 0;         // simulation time
double dt = 0.00001;  // time step
bool simulating = true;

// Global values also accessible by the functions in src/*
int dim = 8;
int dim3 = std::pow(dim, 3.0);
double domain = dim;
bool show_v_field = true;

void draw_vector_field() {
    for (int k = 0; k < dim; k++) {  // Put k on outside to optimize memory access of V_field
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                int flat = flat_index(i, j, k);

                Eigen::RowVector3d V_magnitude_vector(V_field_x(flat), V_field_y(flat), V_field_z(flat));
                V_magnitude_vector = V_magnitude_vector / V_magnitude_vector.norm();
                V_magnitude_vector *= 0.5;

                Eigen::RowVector3d V1((domain / (double)dim) * (double)i, (domain / (double)dim) * (double)j, (domain / (double)dim) * (double)k);
                Eigen::RowVector3d V2 = V1 + V_magnitude_vector;

                int half_dim = dim / 2;
                Eigen::RowVector3d offset;
                offset << half_dim, half_dim, half_dim;

                viz.data().add_edges(V1 - offset, V2 - offset, Eigen::RowVector3d(0, 0, 0));
            }
        }
    }
}

bool simulation_callback() {
    while (simulating) {
        // std::cout << "----------------------- ITER -----------------------" << std::endl;
        // P =  Eigen::MatrixXd::Random(100000,3);

        /******** 1. Apply forces ********/
        Eigen::VectorXd V_field_x_1;
        Eigen::VectorXd V_field_y_1;
        Eigen::VectorXd V_field_z_1;

        Eigen::VectorXd f_x(V_field_x.size());
        Eigen::VectorXd f_y(V_field_y.size());
        Eigen::VectorXd f_z(V_field_z.size());

        for (int i = 0; i < dim3; i++) {
            f_x(i) = 0;
            f_y(i) = 0;
            f_z(i) = 0;
        }

        // add an upwards force of 1 to all the cells on the x,y plane
        for (int i = 1; i < dim-1; i++) {
            for (int j = 1; j < dim-1; j++) {
                f_x(flat_index(i, j, 1)) = 1;
                f_y(flat_index(i, j, 1)) = 1;
                f_z(flat_index(i, j, 1)) = 1;
            }
        }

        add_force(V_field_x_1, V_field_x, f_x, dt);
        add_force(V_field_y_1, V_field_y, f_y, dt);
        add_force(V_field_z_1, V_field_z, f_z, dt);

        apply_fixed_boundary_constraint(V_field_x_1, V_field_y_1, V_field_z_1);

        /******** 2. Advect ********/
        Eigen::VectorXd V_field_x_2;
        Eigen::VectorXd V_field_y_2;
        Eigen::VectorXd V_field_z_2;

        advect(
            V_field_x_2, V_field_y_2, V_field_z_2,  // Output vector field
            V_field_x_1, V_field_y_1, V_field_z_1,  // Input vector field
            dt
        );

        /******** 3. Diffuse ********/
        Eigen::VectorXd V_field_x_3;
        Eigen::VectorXd V_field_y_3;
        Eigen::VectorXd V_field_z_3;

        diffuse( 
            V_field_x_3, V_field_y_3, V_field_z_3,  // Output vector field
            V_field_x_2, V_field_y_2, V_field_z_2,  // Input vector field
            dt
        );


        /******* Update Velocity fields with new values ********/
        // V_field_x = V_field_x_2;
        // V_field_y = V_field_y_2;
        // V_field_z = V_field_z_2;

        V_field_x = V_field_x_3;
        V_field_y = V_field_y_3;
        V_field_z = V_field_z_3;

        t += dt;
    }

    return false;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer) {
    viewer.data().clear();

    // Redraw the points on the visualization
    // viewer.data().clear();
    // viewer.data().point_size = 4;
    // viewer.data().set_points(P,Eigen::RowVector3d(0,0,0));

    // Redraw the velocity field on the visualization

    if (show_v_field)
        draw_vector_field();

    return false;
}

int main(int argc, char **argv) {
    viz.core().background_color.setOnes();

    /******** Create the draw callback ********/
    viz.callback_post_draw = &draw_callback;

    /******** Add particles to the visualization ********/
    P = Eigen::MatrixXd::Random(100000, 3);
    // viz.core().is_animating = true;
    // voz.data().clear();
    // viz.data().point_size = 4;
    // v.data().set_points(P,Eigen::RowVector3d(0,0,0));

    /******** Create the initial velocity field and add them to the visualization ********/
    // Each vector in the velocity field is represented by an entry in each of these vectors (one entry per dimension)
    V_field_x = Eigen::VectorXd(std::pow(dim, 3.0));
    V_field_y = Eigen::VectorXd(std::pow(dim, 3.0));
    V_field_z = Eigen::VectorXd(std::pow(dim, 3.0));

    // This creates a velocity field that initially has vectors moving out from some central point
    // (like an explosion from a central point)
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                int flat = flat_index(i, j, k);

                double x = -1 * (dim / 2) + i;
                double y = -1 * (dim / 2) + j;
                double z = -1 * (dim / 2) + k;
                double denom = std::pow(x, 2.0) + std::pow(y, 2.0) + std::pow(z, 2.0);
                if (denom != 0) {
                    V_field_x(flat) = x / denom;
                    V_field_y(flat) = y / denom;
                    V_field_z(flat) = z / denom;
                } else {
                    V_field_x(flat) = 0;
                    V_field_y(flat) = 0;
                    V_field_z(flat) = 0;
                }
            }
        }
    }

    // for (int i = 1; i < dim-1; i++) {
    //     for (int j = 1; j < dim-1; j++) {
    //         V_field_x(flat_index(i, j, 1)) = 0;
    //         V_field_y(flat_index(i, j, 1)) = 0;
    //         V_field_z(flat_index(i, j, 1)) = 1;
    //     }
    // }

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

    viz.launch();
    return 1;
}
