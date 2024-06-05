#ifndef HH_JACOBISOLVER_HH
#define HH_JACOBISOLVER_HH

#include "ParameterHandler.hpp"
#include "LocalSolver.hpp"
#include "InitializeProblem.hpp"
#include "SafeMPI.hpp"

class JacobiSolver
{
    public:
        JacobiSolver(const std::vector<double> M_, const param::ParameterHandler &params_, int argc_, char **argv_)
        : 
        M(M_), params(params_), argc(argc_), argv(argv_)
        {};

        std::tuple<std::vector<double>, double, unsigned int> solve()
        {
            MPI_Init(&argc, &argv);

            int rank, size;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            // define some parameters just for convenience
            const int N = params.coefficients.n;
            // const double h = 1.0 / N;
            unsigned int n_iterations = 0;
            const unsigned int max_it = params.coefficients.max_it;
            const double tol = params.coefficients.tol_res;
            double norm = std::numeric_limits<double>::infinity();
            // initialize the matrix
            initialize_problem(params, M);

            for(std::size_t k = 0; k < max_it && norm > tol; ++k)
            {
            
                // define the right dimensions for the matrices
                int local_n_cols = N;
                int local_n_rows = (rank < N % size) ? N / size + 1 : N / size; // number of rows for each process
                std::vector<int> local_to_global(size, 0); // local to global index for the first row

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
                auto [local_solution, local_norm] = local_solver.solve();

                // to send data remove the ghost rows
                if(rank == 0)
                {
                    local_solution.erase((local_solution.end() - 1) - local_n_cols, local_solution.end() - 1);
                }
                else if(rank == size - 1)
                {
                    local_solution.erase(local_solution.begin(), local_solution.begin() + local_n_cols);
                }
                else
                {
                    local_solution.erase(local_solution.begin(), local_solution.begin() + local_n_cols);
                    local_solution.erase((local_solution.end() - 1) - local_n_cols, local_solution.end() - 1);
                }

                MPI_Barrier(MPI_COMM_WORLD);

                std::vector<int> send_counts(size, 0);
                std::vector<int> send_start_idx(size, 0);
                std::vector<int> recv_counts(size, 0);
                std::vector<int> recv_start_idx(size, 0);
                
                int start_idx = 0;
                for (int i = 0; i < size; ++i)
                {
                    recv_counts[i] = (N % size > i) ? N / size + 1 : N / size;
                    send_counts[i] = recv_counts[i] * N;

                    recv_start_idx[i] = start_idx;
                    send_start_idx[i] = start_idx * N;

                    start_idx += recv_counts[i];
                }

                std::vector<double> result(N * N, 0.0);

                MPI_Gatherv(local_solution.data(), local_solution.size(), MPI_DOUBLE, 
                            result.data(),  recv_counts.data(), 
                            recv_start_idx.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
                M = result;

                // obtain the norm
                MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                ++n_iterations;
            }


            MPI_Finalize();

            return std::make_tuple(M, norm, n_iterations);
        }

    private:
        std::vector<double> M; // matrix to solve, notice that is just a copy
        param::ParameterHandler params; // parameters
        int argc; // parameter for MPI
        char **argv; // parameter for MPI
    
};

#endif // HH_JACOBISOLVER_HH