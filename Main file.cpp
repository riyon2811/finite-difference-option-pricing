#include "../include/BlackScholesFDM.hpp"
#include <iostream>
#include <iomanip>

/*
    ========================================================================
    Main Driver File
    ------------------------------------------------------------------------

    This file performs the numerical experiment for:

        (a) European Put pricing
        (b) Stability & convergence validation (conceptual discussion)
        (c) Down-and-Out barrier put pricing

    The solver uses the Implicit Euler (Backward Time, Central Space)
    finite difference scheme to approximate the solution of the
    Black–Scholes PDE under the risk-neutral measure.

    Stability Properties of the Scheme:
        - Unconditionally stable (no CFL condition required)
        - First-order accurate in time (O(dt))
        - Second-order accurate in space (O(dS²))
        - Convergent as dt → 0 and dS → 0

    ========================================================================
*/

int main() {

    /*
        --------------------------------------------------------------------
        Model Parameters (Financial Inputs)
        --------------------------------------------------------------------

        S0     : Current stock price
        K      : Strike price
        T      : Time to maturity (in years)
        r      : Constant risk-free interest rate
        sigma  : Volatility of the underlying asset

        These define the Black–Scholes model under the
        risk-neutral probability measure.
    */

    double S0 = 40.0;      // Initial stock price
    double K = 40.0;       // Strike price (ATM option)
    double T = 1.0;        // 1 year to maturity
    double r = 0.10;       // 10% risk-free rate
    double sigma = 0.30;   // 30% annual volatility

    /*
        --------------------------------------------------------------------
        Numerical Grid Parameters
        --------------------------------------------------------------------

        S_max : Upper truncation boundary for stock price domain.
                Since S ∈ [0, ∞), we approximate infinity by S_max.

        M     : Number of spatial grid steps.
                Controls stock price resolution (ΔS).

        N     : Number of time steps.
                Controls temporal resolution (Δt).

        Finer grids (larger M, N) reduce discretization error
        and improve convergence toward analytical solution.
    */

    double S_max = 100.0;  // Upper boundary for stock price grid
    int M = 400;           // Spatial steps
    int N = 400;           // Time steps

    /*
        --------------------------------------------------------------------
        Barrier Parameter (for Down-and-Out Option)
        --------------------------------------------------------------------

        barrier : Knock-out level B.
                  If stock price falls to or below B,
                  the option immediately becomes worthless.

        In this experiment:
            B = 30
    */

    double barrier = 30.0;

    /*
        --------------------------------------------------------------------
        Model Construction
        --------------------------------------------------------------------

        Initializes the finite difference solver with:

            - Model parameters
            - Numerical grid parameters

        Internally computes:
            dS = S_max / M
            dt = T / N
    */

    BlackScholesFDM model(S_max, K, T, r, sigma, M, N);

    /*
        --------------------------------------------------------------------
        Pricing Computations
        --------------------------------------------------------------------

        priceEuropeanPut(S0):
            Solves the Black–Scholes PDE backward in time
            with standard boundary conditions.

        priceDownAndOutPut(S0, barrier):
            Solves modified PDE with absorbing barrier
            condition enforced at S ≤ B.

        Both computations use:
            - Implicit Euler discretization
            - Thomas algorithm for tridiagonal systems
    */

    double europeanPut = model.priceEuropeanPut(S0);
    double barrierPut  = model.priceDownAndOutPut(S0, barrier);

    /*
        --------------------------------------------------------------------
        Output Formatting
        --------------------------------------------------------------------

        std::fixed           : Use fixed-point notation
        std::setprecision(6) : Display 6 decimal places

        This ensures consistent numerical reporting
        for comparison and validation.
    */

    std::cout << std::fixed << std::setprecision(6);

    /*
        --------------------------------------------------------------------
        Display Results
        --------------------------------------------------------------------

        European Put:
            Should converge to analytical Black–Scholes price
            as M and N increase.

        Down-and-Out Put:
            Should be strictly lower than vanilla put price
            due to knock-out risk.

        These values can be used to:
            - Study convergence behavior
            - Analyze discretization bias
            - Compare with closed-form solution
    */

    std::cout << "European Put Price (Implicit Euler FDM): "
              << europeanPut << std::endl;

    std::cout << "Down-and-Out Put Price (Implicit Euler FDM): "
              << barrierPut << std::endl;

    return 0;
}