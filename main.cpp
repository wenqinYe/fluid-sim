#include <iostream>
#include <thread>
#include <visualization.h>
#include <add_force.h>
#include <common.h>
#include <advect.h>
#include <cmath>

igl::opengl::glfw::Viewer viz;

// Simulation state 
Eigen::MatrixXd P;
Eigen::VectorXd V_field_x;
Eigen::VectorXd V_field_y;
Eigen::VectorXd V_field_z;

double t = 0; //simulation time 
double dt = 0.00001; //time step
bool simulating = true;

// Global values also accessible by the functions in src/*
int dim = 10;
int dim3 = std::pow(dim, 3.0);
double domain = 5;

bool simulation_callback() {
    while(simulating) {
        //P =  Eigen::MatrixXd::Random(100000,3);

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
        
        // add a diagonal upswards force of 1 to the first few velocity cells
        for (int i = 2; i < 8; i++) {
            for (int j = 2; j < 8; j++) {
                // f_x(flat_index(i, j, 0)) = 0; 
                // f_y(flat_index(i, j, 0)) = 0; 
                // f_z(flat_index(i, j, 0)) = 1; 
            }
        }

        add_force(V_field_x_1, V_field_x, f_x, dt);
        add_force(V_field_y_1, V_field_y, f_y, dt);
        add_force(V_field_z_1, V_field_z, f_z, dt);


        /******** 2. Advect ********/
        Eigen::VectorXd V_field_x_2;
        Eigen::VectorXd V_field_y_2;
        Eigen::VectorXd V_field_z_2;

        advect(
            V_field_x_2, V_field_y_2, V_field_z_2, // Output vector field
            V_field_x_1, V_field_y_1, V_field_z_1, // Input vector field
            dt
        );

        V_field_x = V_field_x_2;
        V_field_y = V_field_y_2;
        V_field_z = V_field_z_2;

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

    bool show_v_field = true;
    if (show_v_field) {
        for (int k = 0; k < dim; k++) { // Put k on outside to optimize memory access of V_field
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    int flat = flat_index(i, j, k);

                    Eigen::RowVector3d V_magnitude_vector(V_field_x(flat), V_field_y(flat), V_field_z(flat));

                    Eigen::RowVector3d V1((domain/(double)dim)* (double)i, (domain/(double)dim)*(double)j,(domain/(double)dim)*(double)k);
                    Eigen::RowVector3d V2 = V1 + V_magnitude_vector;

                    viz.data().add_edges(V1, V2, Eigen::RowVector3d(0,0,0));
                }
            }
        }
    }

    return false;
}

int main(int argc, char **argv) {
    /******** Create the draw callback ********/
    viz.callback_post_draw = &draw_callback;

    /******** Add particles to the visualization ********/
    P =  Eigen::MatrixXd::Random(100000,3);
    // viz.core().is_animating = true;
    // voz.data().clear();
    // viz.data().point_size = 4;
    // v.data().set_points(P,Eigen::RowVector3d(0,0,0));

    /******** Create the initial velocity field and add them to the visualization ********/
    // Each vector in the velocity field is represented by an entry in each of these vectors (one entry per dimension)
    V_field_x = Eigen::VectorXd(std::pow(dim, 3.0));
    V_field_y = Eigen::VectorXd(std::pow(dim, 3.0));
    V_field_z = Eigen::VectorXd(std::pow(dim, 3.0));

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    int flat = flat_index(i, j, k);
                    V_field_x(flat) = i / (double) dim;
                    V_field_y(flat) = j / (double) dim;
                    V_field_z(flat) = k / (double) dim;
                }
        }
    }

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
    bool show_v_field = true;
    if (show_v_field) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    int flat = flat_index(i, j, k);

                    Eigen::RowVector3d V_magnitude_vector(V_field_x(flat), V_field_y(flat), V_field_z(flat));
                    
                    Eigen::RowVector3d V1((domain/(double)dim)* (double)i, (domain/(double)dim)*(double)j,(domain/(double)dim)*(double)k);
                    Eigen::RowVector3d V2 = V1 + V_magnitude_vector;

                    viz.data().add_edges(V1, V2, Eigen::RowVector3d(0,0,0));
                }
            }
        }
    }

    /******** Start simulation on background thread ********/
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    viz.launch();
    return 1; 

}
