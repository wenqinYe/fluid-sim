#include <add_force.h>
#include <cmath>
#include <math.h>
#include <iostream>

void add_force(Eigen::VectorXd &w1, Eigen::VectorXd &w0, Eigen::VectorXd f, double dt) {
   w1 = w0 + dt * f;
}