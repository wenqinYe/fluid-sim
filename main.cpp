#include <iostream>
#include <thread>
#include <visualization.h>
#include <add_force.h>

igl::opengl::glfw::Viewer viz;
Eigen::MatrixXd P;
double t = 0; //simulation time 
double dt = 0.00001; //time step
int dim = 10;
bool simulating = true;

bool simulation_callback() {

    while(simulating) {
        // Move the particles 
        //P =  Eigen::MatrixXd::Random(100000,3);


    }

    return false;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer) {
    // Redraw the points on the visualization
    // viewer.data().clear();
    // viewer.data().point_size = 4;
    // viewer.data().set_points(P,Eigen::RowVector3d(0,0,0));

    // Redraw the velocity field on the visualization

    return false;
}

int flat_index(int i, int j, int k) {
    return (i * dim) + j + (dim * dim * k); 
}


int main(int argc, char **argv) {
    /******** Start simulation on background thread ********/
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    /******** Create the draw callback ********/
    viz.callback_post_draw = &draw_callback;

    /******** Add particles to the visualization ********/
    P =  Eigen::MatrixXd::Random(100000,3);
    // viz.core().is_animating = true;
    // voz.data().clear();
    // viz.data().point_size = 4;
    // v.data().set_points(P,Eigen::RowVector3d(0,0,0));

    /******** Create the initial velocity field and add them to the visualization ********/
    double domain = 1;

    // Each vector in the velocity field is represented by an entry in each of these vectors (one entry per dimension)
    Eigen::VectorXd V_field_x(std::pow(dim, 3.0));
    Eigen::VectorXd V_field_y(std::pow(dim, 3.0));
    Eigen::VectorXd V_field_z(std::pow(dim, 3.0));

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    int flat = flat_index(i, j, k);
                    V_field_x(flat) = 0.1;
                    V_field_y(flat) = 0.1;
                    V_field_z(flat) = 0.1;
                }
        }
    }

    // Show the velocity_field  in the visualizaiton
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

    viz.launch();
    return 1; 

}
