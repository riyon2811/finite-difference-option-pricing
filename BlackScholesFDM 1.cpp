#ifndef BLACK_SCHOLES_FDM_HPP
#define BLACK_SCHOLES_FDM_HPP

#include <vector>

/*
    ========================================================================
    BlackScholesFDM
    ------------------------------------------------------------------------
    Finite Difference solver for the Black–Scholes PDE under the
    risk-neutral measure.

    This class implements the Implicit Euler (Backward Time, Central Space)
    finite difference scheme for pricing:

        1. European Put Options
        2. Down-and-Out Barrier Put Options

    The resulting linear system at each time step is tridiagonal and
    solved using the Thomas algorithm.

    Numerical Characteristics:
        - Unconditionally stable (Implicit scheme)
        - First-order accurate in time
        - Second-order accurate in space
        - Computational complexity: O(M × N)

    ========================================================================
*/

class BlackScholesFDM {
private:

    /*
        -------------------------
        Model Parameters
        -------------------------

        S_max : Upper truncation boundary for stock price.
                Since S ∈ [0, ∞), we approximate infinity by S_max.

        K     : Strike price of the option.

        T     : Time to maturity (in years).

        r     : Constant risk-free interest rate.

        sigma : Volatility of the underlying asset.

        M     : Number of spatial grid steps (stock price discretization).

        N     : Number of time steps (temporal discretization).

        dS    : Spatial step size = S_max / M.

        dt    : Time step size = T / N.
    */

    double S_max;
    double K;
    double T;
    double r;
    double sigma;
    int M;
    int N;
    double dS;
    double dt;

    /*
        --------------------------------------------------------------------
        solve_tridiagonal(...)
        --------------------------------------------------------------------

        Solves a tridiagonal linear system of the form:

            A x = d

        where A is tridiagonal with:
            a -> sub-diagonal coefficients
            b -> main diagonal coefficients
            c -> super-diagonal coefficients
            d -> right-hand side vector

        This function implements the Thomas algorithm:
            - Forward elimination
            - Backward substitution

        Complexity:
            O(M) per time step

        Purpose in PDE solver:
            At each backward time step, the implicit discretization
            produces a tridiagonal linear system that must be solved
            to obtain the option value vector at the current time level.

        Returns:
            Solution vector x.
    */
    std::vector<double> solve_tridiagonal(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& d
    ) const;

public:

    /*
        --------------------------------------------------------------------
        Constructor
        --------------------------------------------------------------------

        Initializes the Black–Scholes finite difference model.

        Parameters:
            S_max_  : Truncated upper boundary for stock price.
            K_      : Strike price.
            T_      : Time to maturity (years).
            r_      : Risk-free rate.
            sigma_  : Volatility.
            M_      : Number of spatial grid steps.
            N_      : Number of time steps.

        Upon construction:
            dS = S_max / M
            dt = T / N

        These define the discretized computational grid.
    */
    BlackScholesFDM(double S_max_,
                    double K_,
                    double T_,
                    double r_,
                    double sigma_,
                    int M_,
                    int N_);

    /*
        --------------------------------------------------------------------
        priceEuropeanPut(S0)
        --------------------------------------------------------------------

        Prices a European Put Option using the Implicit Euler scheme.

        Steps:
            1. Initialize terminal payoff:
               V(S, T) = max(K - S, 0)

            2. Iterate backward in time:
               Solve A V^n = V^(n+1)

            3. Apply boundary conditions:
               - V(0, t) = K e^{-r(T - t)}
               - V(S_max, t) = 0

            4. Extract option value at S = S0.

        Parameter:
            S0 : Current stock price.

        Returns:
            Numerical approximation of European put price.
    */
    double priceEuropeanPut(double S0) const;

    /*
        --------------------------------------------------------------------
        priceDownAndOutPut(S0, B)
        --------------------------------------------------------------------

        Prices a Down-and-Out Barrier Put Option.

        Additional Feature:
            Absorbing barrier at level B.

        Barrier Condition:
            If S ≤ B:
                Option value = 0
            (Option is knocked out.)

        Implementation Notes:
            - Terminal payoff modified to reflect barrier.
            - At each time step, enforce:
                V(S ≤ B, t) = 0
            - Boundary condition at S=0 becomes absorbing.

        Parameters:
            S0 : Current stock price.
            B  : Barrier level.

        Returns:
            Numerical approximation of barrier put price.
    */
    double priceDownAndOutPut(double S0, double B) const;
};

#endif