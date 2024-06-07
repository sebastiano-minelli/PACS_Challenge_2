#ifndef HH_JACOBISOLVER_HH
#define HH_JACOBISOLVER_HH

#include "ParameterHandler.hpp"
#include "LocalSolver.hpp"
#include "InitializeProblem.hpp"
#include "SafeMPI.hpp"

class JacobiSolver
{
    private:
        std::vector<double> M; // matrix to solve
        param::ParameterHandler params; // parameters
        int argc; // parameter for MPI
        char **argv; // parameter for MPI
        int rank; // MPI rank
        int size; // MPI size

    public:
        JacobiSolver(const std::vector<double> &M_, const param::ParameterHandler &params_, int argc_, char **argv_, int rank_, int size_)
        : 
        M(M_), 
        params(params_), 
        argc(argc_), 
        argv(argv_), 
        rank(rank_), 
        size(size_)
        {};

        // solve function, it returns a tuple with the solution, the norm and the number of iterations
        std::tuple<std::vector<double>, double, unsigned int> 
        solve()
        {
            // define some parameters just for convenience
            const int N = params.coefficients.n; // number of rows (and columns)
            const double h = 1.0 / N; // step size
            unsigned int local_n_iterations = 0; // number of iterations for each process
            const unsigned int max_it = params.coefficients.max_it; // maximum number of iterations
            const double tol = params.coefficients.tol_res; // tolerance for the residual
            double norm = std::numeric_limits<double>::infinity(); // norm of the residual (the global one)
            double local_norm = std::numeric_limits<double>::infinity(); // norm of the residual (the local one)
            std::vector<double> result(N * N, 0.0); // the final solution

            // initialize the matrix (i.e. apply BCs)
            initialize_problem(params, M);

            // define the right dimensions for the local matrix
            int local_n_cols = N;
            int local_n_rows = (rank < N % size) ? N / size + 1 : N / size; // number of rows for each process
            std::vector<int> local_to_global(size, 0); // local to global index, it returns the global matrix first index owned by each process

            // Calculate the number of rows for each process
            int rows_per_process = N / size;
            int remainder = N % size;

            for(int i = 0; i < size; ++i)
            {
                local_to_global[i] = i * rows_per_process + std::min(i, remainder);

                // accounting for the ghost rows (basically every process except the first one has a ghost row on top)
                if(local_to_global[i] != 0)
                    local_to_global[i] -= 1;
            }
            
            // now finish the number of rows calculation
            // add the ghost rows
            if (rank == 0 || rank == size - 1)
            {
                if(local_n_rows < N) // if we run sequentially we don't need to add the ghost rows
                    local_n_rows += 1;
            }
            else // if(rank != 0 && rank != size - 1)
            {
                if(local_n_rows < N - 1 ) // just to make sure to not exceed. Actually for real application size < N almost certainly
                    local_n_rows += 2;
                else
                    local_n_rows += 1;
            } //else do nothing

            
            // define the local matrix
            std::vector<double> local_solution(local_n_rows * local_n_cols, 0.0);

            // Fill the local matrix with the right values
            for(int i = 0; i < local_n_rows; ++i)
            {
                for(int j = 0; j < local_n_cols; ++j)
                {
                    local_solution[i * local_n_cols + j] = M[(local_to_global[rank] + i) * N + j];
                }
            }

            // define the local solver class
            LocalSolver local_solver(local_n_rows, local_n_cols, params, local_to_global[rank]);

            // start the iterations
            for(std::size_t k = 0; k < max_it && local_norm > tol; ++k)
            {
                // update the local solution and get the norm            
                std::tie(local_solution, local_norm) = local_solver.solve(local_solution);

                /*
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
                */
                //std::cout << "End of for loop " << k << " from process " << rank << std::endl;

                if(size > 1)
                {
                    // to gather data remove the ghost rows   
                    if(rank == 0)
                    {
                        // remove the last row
                        local_solution.erase(local_solution.end() - local_n_cols, local_solution.end());
                    }
                    else if(rank != 0 && rank == size - 1)
                    {
                        // remove the first row
                        local_solution.erase(local_solution.begin(), local_solution.begin() + local_n_cols);
                    }
                    else // if(rank != 0 && rank != size - 1)
                    {
                        // remove the first and the last row
                        local_solution.erase(local_solution.end() - local_n_cols, local_solution.end());
                        local_solution.erase(local_solution.begin(), local_solution.begin() + local_n_cols);
                    }
                }

                // send the number of elements to every process
                std::vector<int> recv_counts(size, 0);
                int send_size = local_solution.size();

                MPI_Allgather(&send_size, 1, MPI_INT, 
                            recv_counts.data(), 1, MPI_INT, 
                            MPI_COMM_WORLD);

                // compute the cumulative index for the receive vector buffer
                std::vector<int> recv_start_idx(size, 0); // recv_start_idx[0] = 0
                for (int i = 1; i < size; ++i)
                    recv_start_idx[i] = recv_start_idx[i - 1] + recv_counts[i - 1];

                // gather the data (assemble the global solution matrix)
                MPI_Allgatherv(local_solution.data(), send_size, MPI_DOUBLE, 
                            result.data(),  recv_counts.data(), recv_start_idx.data(), MPI_DOUBLE, 
                            MPI_COMM_WORLD);
                
                // Now that we have the global solution matrix we can subdivide it again and update the local matrix
                
                // Calculate the number of rows for each process
                // We have already computed the local number of rows and columns, no need to do it again (same for the local_to_global vector)

                // Resize the local matrix (that was shrunk before)
                local_solution.clear();
                local_solution.resize(local_n_rows * local_n_cols, 0.0);

                // Subdivide the matrix along different processes
                for(int i = 0; i < local_n_rows; ++i)
                {
                    for(int j = 0; j < local_n_cols; ++j)
                    {
                        local_solution[i * local_n_cols + j] = result[(local_to_global[rank] + i) * N + j];
                    }
                }

                // update local iterations
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

            // notice that we already computed the correct result vector inside the loop
            return std::make_tuple(result, norm, n_iterations);
        }    
};

#endif // HH_JACOBISOLVER_HH