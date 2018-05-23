#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0] = 0.0;

    for (int i = 0; i < N-1; i++) {
      fg[0] += w_cte * CppAD::pow(vars[cte_start + i], 2);
      fg[0] += w_epsi * CppAD::pow(vars[epsi_start + i], 2);
      fg[0] += w_v * CppAD::pow(vars[v_start + i] - max_velocity, 2);
    }

    // Minimize the use of actuators.
    for (int i = 0; i < N - 1 ; i++) {
      fg[0] += w_delta * CppAD::pow(vars[delta_start + i], 2);
      fg[0] += w_a * CppAD::pow(vars[a_start + i], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int i = 0; i < N - 2; i++) {
      fg[0] += w_ddelta * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += w_da * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    // Initial constraints given state does not vary
    fg[x_start + 1] = vars[x_start];
    fg[y_start + 1] = vars[y_start];
    fg[psi_start + 1] = vars[psi_start];
    fg[v_start + 1] = vars[v_start];
    fg[cte_start + 1] = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];

    // constraints based on our kinematic model
    for (int i = 0; i < N - 1; ++i) {
      //current state and actuations
      const auto px0 = vars[x_start + i];
      const auto py0 = vars[y_start + i];
      const auto psi0 = vars[psi_start + i];
      const auto v0 = vars[v_start + i];
      const auto cte0 = vars[cte_start + i];
      const auto epsi0 = vars[epsi_start + i];
      const auto delta0 = vars[delta_start + i];
      const auto a0 = vars[a_start + i];

      // next state
      const auto px1 = vars[x_start + i + 1];
      const auto py1 = vars[y_start + i + 1];
      const auto psi1 = vars[psi_start + i + 1];
      const auto v1 = vars[v_start + i + 1];
      const auto cte1 = vars[cte_start + i + 1];
      const auto epsi1 = vars[epsi_start + i + 1];

      // desired py and psi
      const auto py_desired = coeffs[3] * CppAD::pow(px0, 3) + coeffs[2] * CppAD::pow(px0, 2) + coeffs[1] * px0 + coeffs[0];
      const auto psi_desired = CppAD::atan(3.0 * coeffs[3] * CppAD::pow(px0, 2) + 2.0 * coeffs[2] * px0 + coeffs[1]);

      // relationship of current state + actuations and next state
      // based on our kinematic model
      const auto px1_f = px0 + v0 * CppAD::cos(psi0) * dt;
      const auto py1_f = py0 + v0 * CppAD::sin(psi0) * dt;
      const auto psi1_f = psi0 + v0 * (-delta0) / Lf * dt;
      const auto v1_f = v0 + a0 * dt;
      const auto cte1_f = py_desired - py0 + v0 * CppAD::sin(epsi0) * dt;
      const auto epsi1_f = psi0 - psi_desired + v0 * (-delta0) / Lf * dt;

      // store the constraint expression of two consecutive states
      fg[x_start + i + 2] = px1 - px1_f;
      fg[y_start + i + 2] = py1 - py1_f;
      fg[psi_start + i + 2] = psi1 - psi1_f;
      fg[v_start + i + 2] = v1 - v1_f;
      fg[cte_start + i + 2] = cte1 - cte1_f;
      fg[epsi_start + i + 2] = epsi1 - epsi1_f;
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {
  this->vars.resize(number_of_states_actuations);
  this->vars_lowerbound.resize(number_of_states_actuations);
  this->vars_upperbound.resize(number_of_states_actuations);
  this->constraints_lowerbound.resize(number_of_constraints);
  this->constraints_upperbound.resize(number_of_constraints);

  for (int i = 0; i < number_of_states_actuations; i++) {
    this->vars[i] = 0.0;
  }

  for (int i = 0; i < number_of_constraints; i++) {
    this->constraints_lowerbound[i] = 0.0;
    this->constraints_upperbound[i] = 0.0;
  }
  
  for (int i = 0; i < delta_start; i++) {
    this->vars_lowerbound[i] = -1.0e10;
    this->vars_upperbound[i] = 1.0e10;
  }

  // Actuation variables i.e. steering acceleration should be between [-1, 1]
  for (int i = delta_start; i < a_start; i++) {
    this->vars_lowerbound[i] = -0.75;
    this->vars_upperbound[i] = 0.75;
  }

  for (int i = a_start; i < number_of_states_actuations; i++) {
    this->vars_lowerbound[i] = -0.5;
    this->vars_upperbound[i] = 1.0;
  }

}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:

  this->vars[x_start] = state[0];
  this->vars[y_start] = state[1];
  this->vars[psi_start] = state[2];
  this->vars[v_start] = state[3];
  this->vars[cte_start] = state[4];
  this->vars[epsi_start] = state[5];

  this->constraints_lowerbound[x_start] = state[0];
  this->constraints_lowerbound[y_start] = state[1];
  this->constraints_lowerbound[psi_start] = state[2];
  this->constraints_lowerbound[v_start] = state[3];
  this->constraints_lowerbound[cte_start] = state[4];
  this->constraints_lowerbound[epsi_start] = state[5];

  this->constraints_upperbound[x_start] = state[0];
  this->constraints_upperbound[y_start] = state[1];
  this->constraints_upperbound[psi_start] = state[2];
  this->constraints_upperbound[v_start] = state[3];
  this->constraints_upperbound[cte_start] = state[4];
  this->constraints_upperbound[epsi_start] = state[5];

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  this->steer = solution.x[delta_start];
  this->throttle = solution.x[a_start];
  this->future_x = {};
  this->future_y = {};
  for (int i = 0; i < N; i++) {
    this->future_x.emplace_back(solution.x[x_start + i]);
    this->future_y.emplace_back(solution.x[y_start + i]);
  }
  vector<double> result;
  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);
  return result;
}
