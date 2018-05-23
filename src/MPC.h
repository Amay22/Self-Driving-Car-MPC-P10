#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <cppad/cppad.hpp>

typedef CPPAD_TESTVECTOR(double) Dvector;

using namespace std;

// timestep length and duration
const int N = 10;
const double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

const double max_velocity = 100.0;
const int number_of_states = 6; // px, py, psi, v, cte, epsi
const int number_of_states_actuations =  N * number_of_states + (N - 1) * 2;
const int number_of_constraints = N * number_of_states;

const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

// weights for cost computations
const double w_cte = 1500.0;
const double w_epsi = 1500.0;
const double w_v = 1.0;
const double w_delta = 10.0;
const double w_a = 10.0;
const double w_ddelta = 150.0;
const double w_da = 15.0;

class MPC {
 public:

  double steer;
  double throttle;

  Dvector vars; // where all the state and actuation variables will be stored
  Dvector vars_lowerbound; //lower limit for each corresponding variable in x
  Dvector vars_upperbound; //upper limit for each corresponding variable in x
  Dvector constraints_lowerbound; // value constraint for each corresponding constraint expression
  Dvector constraints_upperbound; // value constraint for each corresponding constraint expression

  std::vector<double> future_x;
  std::vector<double> future_y;

  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
