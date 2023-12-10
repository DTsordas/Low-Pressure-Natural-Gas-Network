/// Included libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

/// Definition of all MACROS
#define N_I 1.15
#define Roughness 0.0000115
#define Threshold 0.0001
#define Kinematic_Viscosity 0.00001357
#define Density 0.79
#define MAX_ITER 10000

/// Definition of sign Function
int sign(double x) {
    if(x >= 0) {
        return 1;
    } else {
        return -1;
    }
}

/// Definition of LU Solver
void LU(const std::vector<std::vector<double>>& A, std::vector<double>& X, const std::vector<double>& B) {
    int n = static_cast<int>(A.size());
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));

    /// Doolittle's LU decomposition
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double sum = 0.0;

            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }

            U[i][j] = A[i][j] - sum;
        }

        for (int j = i; j < n; j++) {
            if (i == j) {
                L[i][j] = 1.0;
            } else {
                double sum = 0.0;

                for (int k = 0; k < i; k++) {
                    sum += L[j][k] * U[k][i];
                }

                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }
    }

    /// Forward substitution (Ly = B)
    std::vector<double> Y(n, 0.0);
    for (int i = 0; i < n; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
    }

    /// Backward substitution (Ux = y)
    X.resize(n);
    for (int i = n - 1; i >= 0; i--) {
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
}

int main() {
    /// This map represents the network connections, mapping each node to its connected nodes via each pipe.
    std::vector<std::vector<int>> Nodes{{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {1, 0, 2, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 2, 0, 3, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 3, 0, 4, 0, 0, 7, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 5, 0, 0, 0, 0, 8, 0, 10, 0, 0, 0, 0, 0},
                                        {0, 0, 6, 0, 0, 8, 0, 9, 0, 11, 0, 0, 0, 0},
                                        {0, 0, 0, 7, 0, 0, 9, 0, 0, 0, 12, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 10, 0, 0, 0, 13, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 11, 0, 13, 0, 14, 15, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 12, 0, 14, 0, 0, 16, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 17, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 17, 0, 18},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0}
    };

    int const Pipes = 18;

    /// Vectors to store the initial data.
    std::vector<double> Pnodes(Nodes.size(), 0.0);       /// Pressure (Pa) -> [14]
    std::vector<double> Pstatic(Nodes.size(), 0.0);      /// Static Pressure (Pa) -> [14]
    std::vector<double> QLoses(Nodes.size(), 0.0);       /// Loses in each Node (m^3/h) -> [14]
    std::vector<double> k(Pipes,0.0);                    /// Linear loses -> [18]
    std::vector<double> Lamda(Pipes,0.0);                /// Friction Factor -> [18]
    std::vector<double> Qpipe(Pipes,0.0);                /// Flow Rate (m^3/h) -> [18]
    std::vector<double> dpPipes(Pipes,0.0);              /// Total Δh(Pa) -> [18]
    std::vector<double> Length(Pipes,0.0);               /// Length (m) -> [18]
    std::vector<double> Diameter(Pipes,0.0);             /// Diameter (m) -> [18]
    std::vector<double> Re(Pipes,0.0);                   /// Reynolds number -> [18]
    std::vector<double> Velocity(Pipes,0.0);             /// Velocity of each pipe (m/sec) -> [18]
    std::vector<double> Z(Pipes,0.0);                    /// ζ (sec^2/m^5) -> [18]

    /// Input Data.
    {
        Pnodes = {3000.0, 2980.0, 2930.0, 2900.0, 2880.0, 2810.0, 2760.0, 2740.0, 2700.0, 2680.0, 2650.0, 2620.0 ,2580.0, 2550.0};
        QLoses = {1000, -40, -60, -80, -150, -70, -70, -60, -20, -40, -50, -80, -60, -220};
        Length = {60, 60, 60, 40, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 40};
        /// Initial Diameters
        //  Diameter = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};
        /// Optimized Diameters
        //  Diameter = {0.363575 ,0.253739 ,0.181075 ,0.140831 ,0.250032 ,0.153836 ,0.0488543 ,0.165815 ,0.144144 ,0.160521 ,0.145354 ,0.123406 ,0.152062 ,0.121563 ,0.155517 ,0.152972 ,0.116663 ,0.170531};
        /// Optimized SDR Diameters
        Diameter = {0.3682 ,0.2292 ,0.1840 ,0.1472 ,0.2578 ,0.1636 ,0.0514 ,0.1840 ,0.1472 ,0.1636 ,0.1472 ,0.1308 ,0.1636 ,0.1308 ,0.1636 ,0.1636 ,0.1146 ,0.1840};
        Z = {0 ,0.7 ,0.7 ,0.7 ,1.7 ,1.7 ,1.7 ,1.7 ,1.7 ,0.7 ,0.7 ,0.7 ,1.3 ,0.7 ,1.7 ,0.7 ,1.3 ,0.7};
    }

    /// Conversion of QLoses from m^3/h to m^3/sec.
    for (double & i : QLoses) {
        i = (i * N_I * 1.3) / 3600;
    }

    /// Initialize the F vector to store flow corrections for each node.
    std::vector<double> F(Nodes.size(), 0.0);
    /// Initialize the dh vector to store head loss corrections for each node.
    std::vector<double> dh(Nodes.size(), 0.0);
    /// Initialize the Jacobian matrix to be used in the iterative solver.
    std::vector<std::vector<double>> Jacobian(Nodes.size(), std::vector<double>(Nodes.size(), 0.0));

    /// Pre-calculate the friction factor (λ) and the total frictional loss coefficient (k) for each pipe.
    for (int i = 0; i < Pipes; ++i) {
        Lamda[i] = pow(1 / (1.14 - 2 * log10((Roughness / Diameter[i]))), 2);
        k[i] = ((Lamda[i] * Length[i] * 8 * Density) / (M_PI * M_PI * pow(Diameter[i], 5)));
    }

    for (int x = 0; x < MAX_ITER; ++x) {

        /// Compute the Jacobian matrix and F vector for the current iteration.
        for (int i = 0; i < Nodes.size(); ++i) {
            for (int j = 0; j < Nodes.size(); ++j) {
                /// Define the index of the pipe connecting nodes i and j.
                int pipeIndex = Nodes[i][j];
                if (pipeIndex != 0) {
                    /// Update the Jacobian matrix and F vector based on the pressure differences and flow characteristics.
                    dpPipes[pipeIndex - 1] = Pnodes[i] - Pnodes[j];
                    Jacobian[i][0] = 0;
                    Jacobian[i][j]  = (0.5 * (1 / sqrt(k[pipeIndex - 1])) * (1 / sqrt(fabs(dpPipes[pipeIndex - 1]))));
                    Jacobian[i][i] -= (0.5 * (1 / sqrt(k[pipeIndex - 1])) * (1 / sqrt(fabs(dpPipes[pipeIndex - 1]))));
                    F[i] += sign((dpPipes[pipeIndex - 1])) * sqrt(fabs(dpPipes[pipeIndex - 1]) / k[pipeIndex - 1]);
                }
            }
            /// Incorporate the effect of localized losses for each node into the F vector.
            F[i] -= QLoses[i];
        }

        /// Apply boundary conditions specific to the problem being solved.
        Jacobian[0][0] = 1;

        /// Solve the linear system Jacobian * Δh = F using LU decomposition.
        LU(Jacobian, dh, F);

        /// Break Condition
        {
            double RMS = 0;
            for (int i = 0; i < Nodes.size(); ++i) {
                RMS += dh[i] * dh[i];
            }
            RMS = sqrt(RMS) / static_cast<double>(Nodes.size());

            if (RMS <= Threshold) {
                std::cout << "SYSTEM CONVERGED IN " << x << " ITERATIONS." << std::endl;
                std::cout << std::endl;
                break;
            }
        }

        for (int i = 0; i < Nodes.size(); ++i) {
            /// Apply correction to the pressure of each node.
            Pnodes[i]+= dh[i];
            for (int j = 0; j < Nodes.size(); ++j) {
                /// Define the index of the pipe connecting nodes i and j.
                int pipeIndex = Nodes[i][j];
                if (pipeIndex != 0) {
                    /// Calculate the new pressure difference (Δh) between nodes for each pipe.
                    dpPipes[pipeIndex - 1] = Pnodes[i] - Pnodes[j];
                    /// Calculate the flow rate in each pipe.
                    Qpipe[pipeIndex - 1] = sqrt(fabs(dpPipes[pipeIndex - 1]) / k[pipeIndex - 1]);
                    /// Calculate the Reynolds number for each pipe.
                    Re[pipeIndex - 1] = ((4 * Qpipe[pipeIndex - 1]) / (M_PI * Diameter[pipeIndex - 1] * Kinematic_Viscosity));
                    /// Calculate the velocity of flow in each pipe.
                    Velocity[pipeIndex - 1] = ((Re[pipeIndex - 1] * Kinematic_Viscosity) / Diameter[pipeIndex - 1]);
                    /// Calculate the λ friction factor for each pipe.
                    Lamda[pipeIndex - 1] = pow(1 / (1.14 - 2 * log10((Roughness / Diameter[pipeIndex - 1]) + (21.25 / pow(Re[pipeIndex - 1], 0.9)))), 2);
                    /// Calculate the k frictional loss coefficient for each pipe.
                    k[pipeIndex - 1] = (((Lamda[pipeIndex - 1] * Length[pipeIndex - 1]) / Diameter[pipeIndex - 1]) + Z[pipeIndex - 1]) * ((8 * Density) / (M_PI * M_PI * pow(Diameter[pipeIndex - 1], 4)));
                    /// Calculate the Static Pressure for each node.
                    Pstatic[i] = (Pnodes[i] - (Density * Velocity[pipeIndex - 1] * Velocity[pipeIndex - 1] * 0.5));
                }
                /// Initialize the Jacobian matrix elements to zero.
                Jacobian[i][j] = 0;
            }
            /// Reset the dh correction and F vectors for the next iteration.
            F[i] = 0;
            dh[i] = 0;
        }

    }

    const int columnWidth = 25;

    /// Print the header for Total Pressure and Static Pressure.
    std::cout << std::left << std::setw(columnWidth) << "Total Pressure [Pa]" << "Static Pressure [Pa]" << '\n';

    /// Print the values from Pnodes and Pstatic.
    for (int i = 0; i < Nodes.size(); ++i) {
        std::cout << std::left << std::setw(columnWidth) << Pnodes[i]
                  << std::left << std::setw(columnWidth) << Pstatic[i] << '\n';
    }

    std::cout << std::string(75, '-') << '\n';

    /// Print the headers for Velocity, Reynolds, and Q.
    std::cout << std::left << std::setw(columnWidth) << "Velocity [m/sec]"
              << std::left << std::setw(columnWidth) << "Reynolds"
              << std::left << std::setw(columnWidth) << "Q [m^3/sec]" << '\n';

    /// Print the values from Velocity, Re, and Qpipe side by side.
    for (int i = 0; i < Pipes; ++i) { // Since each vector has 18 elements
        std::cout << std::left << std::setw(columnWidth) << Velocity[i]
                  << std::left << std::setw(columnWidth) << Re[i]
                  << std::left << std::setw(columnWidth) << Qpipe[i] << '\n';
    }

    return 0;
}

