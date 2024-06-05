#ifndef HH_JACOBISOLVER_HH
#define HH_JACOBISOLVER_HH

#include "ParameterHandler.hpp"
#include "LocalSolver.hpp"
#include "InitializeProblem.hpp"
#include "SafeMPI.hpp"

class JacobiSolver
{
    public:
        JacobiSolver(const param::ParameterHandler &params_, int argc_, char **argv_)
        : 
        M(M_), params(params_), argc(argc_), argv(argv_)
        {};

        std::tuple<Eigen::MatrixXd, double, unsigned int> solve();
        {
            MPI_Init(&argc, &argv);

            int rank, size;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            // define some parameters just for convenience
            const int N = params.coefficients.n;
            const double h = 1.0 / N;
            const unsigned int max_it = params.solver.max_it;
            const double tol = params.solver.tol;
            
            // define the right dimensions for the matrices
            const int local_n_cols = N;
            const int local_n_rows = (rank < N % size) ? N / size + 1 : N / size; // number of rows for each process
            const std::vector<int> local_to_global(size, 0); // local to global index for the first row

            for(int i = 0; i < size; ++i)
            {
                local_to_global[i] = (i < N % size) ? i * (N / size + 1) : i * (N / size) + N % size;
                // account for the ghost rows
                if(local_to_global[i] != 0)
                    local_to_global[i] -= 1;
            }

            // add the ghost rows
            if (rank == 0 || rank == size - 1)
                local_n_rows += 1;
            else
                local_n_rows += 2;
            
            // define the local matrix
            std::vector<double> M_local(local_n_rows * local_n_cols, 0.0);

            // initialize the matrix
            initialize_problem(params, M);

            // Subdivide the matrix along different processes
            for(int i = 0; i < local_n_rows; ++i)
            {
                for(int j = 0; j < local_n_cols; ++j)
                {
                    M_local[i * local_n_cols + j] = M[(local_to_global[rank] + i) * local_n_cols + j];
                }
            }

            // define the local solution
            LocalSolver local_solver(M_local, params);
            auto [local_solution, norm] = local_solver.solve();

            /*
            std::vector<int> send_counts(size, 0);
            std::vector<int> send_start_idx(size, 0);
            std::vector<int> recv_counts(size, 0);
            std::vector<int> recv_start_idx(size, 0);
            */


            MPI_Finalize();
        }

    private:
        Eigen::MatrixXd M; // matrix to solve
        param::ParameterHandler params; // parameters
        int argc; // parameter for MPI
        char **argv; // parameter for MPI
    
};

#endif // HH_JACOBISOLVER_HH