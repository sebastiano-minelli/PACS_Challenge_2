#ifndef HH_JACOBISOLVER_HH
#define HH_JACOBISOLVER_HH

#include "ParameterHandler.hpp"
#include "LocalSolver.hpp"
#include "InitializeProblem.hpp"
#include "SafeMPI.hpp"

class JacobiSolver
{
    public:
        JacobiSolver(const std::vector<double> M_, const param::ParameterHandler &params_, int argc_, char **argv_, int rank_, int size_)
        : 
        M(M_), params(params_), argc(argc_), argv(argv_), rank(rank_), size(size_)
        {};

        std::tuple<std::vector<double>, double, unsigned int> solve()
        {
            // define some parameters just for convenience
            const int N = params.coefficients.n;
            const double h = 1.0 / N;
            unsigned int local_n_iterations = 0;
            const unsigned int max_it = params.coefficients.max_it;
            const double tol = params.coefficients.tol_res;
            double norm = std::numeric_limits<double>::infinity();

            double local_norm = std::numeric_limits<double>::infinity();
            std::vector<double> result(N * N, 0.0); // the final solution

            // initialize the matrix
            initialize_problem(params, M);

            // define the right dimensions for the matrices
            int local_n_cols = N;
            int local_n_rows = (rank < N % size) ? N / size + 1 : N / size; // number of rows for each process
            std::vector<int> local_to_global(size, 0); // local to global index for the first row

            // Calculate the number of rows for each process
            int rows_per_process = N / size;
            int remainder = N % size;

            for(int i = 0; i < size; ++i)
            {
                local_to_global[i] = i * rows_per_process + std::min(i, remainder);
                // account for the ghost rows
                if(local_to_global[i] != 0)
                    local_to_global[i] -= 1;
            }
            
            // add the ghost rows
            if (rank == 0 || rank == size - 1)
            {
                if(local_n_rows < N)
                    local_n_rows += 1;
            }
            else // if(rank != 0 && rank != size - 1)
            {
                if(local_n_rows < N - 1 )
                    local_n_rows += 2;
                else
                    local_n_rows += 1;
            } //else do nothing

            
            // define the local matrix
            std::vector<double> local_solution(local_n_rows * local_n_cols, 0.0);

            // Subdivide the matrix along different processes
            for(int i = 0; i < local_n_rows; ++i)
            {
                for(int j = 0; j < local_n_cols; ++j)
                {
                    local_solution[i * local_n_cols + j] = M[(local_to_global[rank] + i) * local_n_cols + j];
                }
            }
            // define the local solution
            LocalSolver local_solver(local_n_rows, local_n_cols, params, local_to_global[rank]);

            for(std::size_t k = 0; k < max_it && local_norm > tol; ++k)
            {
            
                std::tie(local_solution, local_norm) = local_solver.solve(local_solution); // local solution is passed by reference and updated

                // passing info to nearby processes
                enum class Tag : int
                {
                    TOP = 0,
                    BOTTOM = 1
                };
                if(size > 1) // if size == 1 we have everything already
                {
                    std::vector<MPI_Request> requests;
                    std::vector<MPI_Status> statuses(2);

                    // initialize receive 
                    if(rank != 0)
                    {
                        requests.emplace_back();
                        // top receive
                        MPI_Irecv(local_solution.data() + N, N, MPI_DOUBLE,
                                rank - 1, static_cast<int>(Tag::BOTTOM), MPI_COMM_WORLD, &requests.back());
                    }

                    // initialize send
                    if(rank != size - 1)
                    {
                        requests.emplace_back();
                        // bottom send
                        MPI_Isend(local_solution.data() + local_solution.size() - 2 * N, N, MPI_DOUBLE,
                                rank + 1, static_cast<int>(Tag::BOTTOM), MPI_COMM_WORLD, &requests.back());
                    }

                    // initialize receive 
                    if(rank != size - 1)
                    {
                        requests.emplace_back();
                        // bottom receive
                        MPI_Irecv(local_solution.data() + local_solution.size() - 2 * N, N, MPI_DOUBLE,
                                rank + 1, static_cast<int>(Tag::BOTTOM), MPI_COMM_WORLD, &requests.back());
                    }

                    // initialize send
                    if(rank != 0)
                    {
                        requests.emplace_back();
                        // top send
                        MPI_Isend(local_solution.data() + N, N, MPI_DOUBLE,
                                rank - 1, static_cast<int>(Tag::BOTTOM), MPI_COMM_WORLD, &requests.back());
                    }

                    // Wait for all requests to complete
                    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
                    // Wait for all non-blocking operations to complete
                    //MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
                    
                }
                //std::cout << "End of for loop " << k << " from process " << rank << std::endl;

                ++local_n_iterations;
            }


            MPI_Barrier(MPI_COMM_WORLD); // wait for all the processes to finish the computation

            // to compute the real norm we need first to recover the sum of the squares of the differences
            local_norm = (local_norm * local_norm) / h;

            // obtain the sum
            MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // compute the total norm
            norm = std::sqrt(h * norm);

            // obtain the max number of iterations
            unsigned int n_iterations;
            MPI_Allreduce(&local_n_iterations, &n_iterations, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

            // to gather data remove the ghost rows   
            if(rank == 0)
            {
                local_solution.erase(local_solution.end() - local_n_cols, local_solution.end());
            }
            else if(rank == size - 1)
            {
                local_solution.erase(local_solution.begin(), local_solution.begin() + local_n_cols);
            }
            else // if(rank != 0 && rank != size - 1)
            {
                local_solution.erase(local_solution.begin(), local_solution.begin() + local_n_cols);
                local_solution.erase(local_solution.end() - local_n_cols, local_solution.end());
            }

            // send the number of elements to every process
            std::vector<int> recv_counts(size, 0);
            int send_size = local_solution.size();
            MPI_Allgather(&send_size, 1, MPI_INT, 
                        recv_counts.data(), 1, MPI_INT, 
                        MPI_COMM_WORLD);

            // compute the displacement
            std::vector<int> recv_start_idx(size, 0); // recv_start_idx[0] = 0
            for (int i = 1; i < size; ++i)
                recv_start_idx[i] = recv_start_idx[i - 1] + recv_counts[i - 1];

            MPI_Allgatherv(local_solution.data(), local_solution.size(), MPI_DOUBLE, 
                        result.data(),  recv_counts.data(), recv_start_idx.data(), MPI_DOUBLE, 
                        MPI_COMM_WORLD);

            return std::make_tuple(result, norm, n_iterations);
        }

    protected:
        std::vector<double> M; // matrix to solve, notice that is just a copy
        param::ParameterHandler params; // parameters
        int argc; // parameter for MPI
        char **argv; // parameter for MPI
        int rank;
        int size;
    
};

#endif // HH_JACOBISOLVER_HH