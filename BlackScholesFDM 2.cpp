#include "../include/BlackScholesFDM.hpp"
#include <cmath>
#include <algorithm>

/*
    ========================================================================
    Constructor
    ------------------------------------------------------------------------
    Initializes model parameters and constructs the computational grid.

    S_max  : Truncated upper boundary for stock price domain.
             Since S ∈ [0, ∞), we approximate infinity by S_max.

    K      : Strike price of the option.

    T      : Time to maturity (in years).

    r      : Constant risk-free interest rate under risk-neutral measure.

    sigma  : Volatility of underlying asset.

    M      : Number of spatial grid steps (discretization of stock price).

    N      : Number of time steps (discretization of time).

    Derived quantities:
        dS = spatial step size
        dt = time step size
    ========================================================================
*/

BlackScholesFDM::BlackScholesFDM(double S_max_,
                                 double K_,
                                 double T_,
                                 double r_,
                                 double sigma_,
                                 int M_,
                                 int N_)
    : S_max(S_max_), K(K_), T(T_),
      r(r_), sigma(sigma_), M(M_), N(N_) {

    dS = S_max / M;   // Spatial grid spacing
    dt = T / N;       // Time step size
}

/*
    ========================================================================
    Thomas Algorithm (Tridiagonal Matrix Algorithm - TDMA)
    ------------------------------------------------------------------------

    Solves linear systems of the form:

        A x = d

    where A is tridiagonal:

        a[i] = sub-diagonal
        b[i] = main diagonal
        c[i] = super-diagonal

    Steps:
        1. Forward elimination
        2. Backward substitution

    Complexity:
        O(M) per solve

    In the context of the PDE:
        At each time step, the implicit discretization produces
        a tridiagonal linear system that must be solved to obtain
        the option value at the previous time level.

    ========================================================================
*/

std::vector<double> BlackScholesFDM::solve_tridiagonal(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d) const {

    int n = b.size();

    std::vector<double> c_star(n);  // Modified super-diagonal
    std::vector<double> d_star(n);  // Modified RHS
    std::vector<double> x(n);       // Solution vector

    // ---- Forward Elimination ----
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {

        // Modify pivot to eliminate sub-diagonal
        double m = b[i] - a[i] * c_star[i - 1];

        c_star[i] = c[i] / m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }

    // ---- Backward Substitution ----
    x[n - 1] = d_star[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    return x;
}

/*
    ========================================================================
    priceEuropeanPut(S0)
    ------------------------------------------------------------------------
    Prices a European Put Option using Implicit Euler (BTCS).

    PDE under risk-neutral measure:

        ∂V/∂t + 0.5 σ² S² ∂²V/∂S² + rS ∂V/∂S − rV = 0

    Terminal condition:
        V(S,T) = max(K − S, 0)

    Boundary conditions:
        V(0,t)     = K e^{−r(T−t)}
        V(S_max,t) = 0

    Algorithm:
        1. Initialize payoff at maturity.
        2. March backward in time.
        3. At each time step, solve tridiagonal system.
        4. Extract value at S = S0.

    ========================================================================
*/

double BlackScholesFDM::priceEuropeanPut(double S0) const {

    // Option value vector across stock price grid
    std::vector<double> V(M + 1);

    // Stock price grid
    std::vector<double> S(M + 1);

    // ---- Terminal Payoff Initialization ----
    for (int i = 0; i <= M; ++i) {
        S[i] = i * dS;
        V[i] = std::max(K - S[i], 0.0);
    }

    // ---- Backward Time Marching ----
    for (int n = N - 1; n >= 0; --n) {

        std::vector<double> a(M + 1), b(M + 1),
                            c(M + 1), d(M + 1);

        // Build coefficient matrix
        for (int i = 1; i < M; ++i) {

            /*
                Discretized coefficients derived from:
                - Central difference in space
                - Backward difference in time
            */

            double alpha = 0.5 * dt *
                (sigma * sigma * i * i - r * i);

            double beta = 1 + dt *
                (sigma * sigma * i * i + r);

            double gamma = -0.5 * dt *
                (sigma * sigma * i * i + r * i);

            a[i] = -alpha;   // Sub-diagonal
            b[i] = beta;     // Main diagonal
            c[i] = -gamma;   // Super-diagonal
            d[i] = V[i];     // RHS from previous time step
        }

        // ---- Boundary Conditions ----

        // At S = 0 (deep in-the-money put)
        b[0] = 1.0;
        d[0] = K * std::exp(-r * (T - n * dt));

        // At S = S_max (deep out-of-the-money put)
        b[M] = 1.0;
        d[M] = 0.0;

        // Solve linear system
        V = solve_tridiagonal(a, b, c, d);
    }

    // Interpolate value at S0 (grid index approximation)
    int idx = static_cast<int>(S0 / dS);
    return V[idx];
}

/*
    ========================================================================
    priceDownAndOutPut(S0, B)
    ------------------------------------------------------------------------
    Prices a Down-and-Out Barrier Put Option.

    Additional Constraint:
        If S ≤ B → option is knocked out (value = 0).

    Implementation:
        - Terminal payoff modified.
        - Enforce absorbing boundary condition at each time step.
        - Barrier introduces discontinuity in solution domain.

    ========================================================================
*/

double BlackScholesFDM::priceDownAndOutPut(double S0, double B) const {

    std::vector<double> V(M + 1);
    std::vector<double> S(M + 1);

    // ---- Terminal Condition with Barrier ----
    for (int i = 0; i <= M; ++i) {

        S[i] = i * dS;

        if (S[i] <= B)
            V[i] = 0.0;  // Knocked out
        else
            V[i] = std::max(K - S[i], 0.0);
    }

    // ---- Backward Time Marching ----
    for (int n = N - 1; n >= 0; --n) {

        std::vector<double> a(M + 1), b(M + 1),
                            c(M + 1), d(M + 1);

        for (int i = 1; i < M; ++i) {

            // Enforce barrier condition
            if (S[i] <= B) {
                b[i] = 1.0;
                d[i] = 0.0;
                continue;
            }

            double alpha = 0.5 * dt *
                (sigma * sigma * i * i - r * i);

            double beta = 1 + dt *
                (sigma * sigma * i * i + r);

            double gamma = -0.5 * dt *
                (sigma * sigma * i * i + r * i);

            a[i] = -alpha;
            b[i] = beta;
            c[i] = -gamma;
            d[i] = V[i];
        }

        // ---- Boundary Conditions ----

        // Lower boundary behaves as absorbing
        b[0] = 1.0;
        d[0] = 0.0;

        // Upper boundary
        b[M] = 1.0;
        d[M] = 0.0;

        V = solve_tridiagonal(a, b, c, d);
    }

    int idx = static_cast<int>(S0 / dS);
    return V[idx];
}