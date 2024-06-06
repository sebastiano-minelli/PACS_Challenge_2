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
            // const double h = 1.0 / N;
            unsigned int n_iterations = 0;
            const unsigned int max_it = params.coefficients.max_it;
            const double tol = params.coefficients.tol_res;
            double norm = std::numeric_limits<double>::infinity();

            double local_norm = 0.0;
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

            for(std::size_t k = 0; k < max_it && norm > tol; ++k)
            {
            
                local_norm = local_solver.solve(local_solution); // local solution is passed by reference and updated
                MPI_Barrier(MPI_COMM_WORLD); // wait for all the processes to finish the computation

                // obtain the norm
                MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // passing info to nearby processes
                /*
                Note about Tags:
                the tag is defined as 1p000ij where:
                p = 1 if we send the last row
                p = 0 if we send the first row
                the receiver has p inverted (0 if the sender sends the last row, 1 if the sender sends the first row)
                i = the process that sends the data
                j = the process that receives the data
                */
                if(size > 1) // if size == 1 we have everything already
                {
                    std::vector<MPI_Request> requests(4*(size-1), MPI_REQUEST_NULL);
                    std::vector<MPI_Status>  statuses(4*(size-1));
                    if(rank == 0)
                    {
                        // send the second last row to process 1 (the last one is the ghost one, process 1 has it already)
                        MPI_Isend(local_solution.data() + local_solution.size() - 2 * N, N, MPI_DOUBLE,
                                1, (1100000 + 0 * 10 + 1), MPI_COMM_WORLD, &requests[0]); 

                        // receive the first row from process 1 (that will be placed in the last ghost row of process 0)
                        MPI_Irecv(local_solution.data() + local_solution.size() - N, N, MPI_DOUBLE,
                                1, (1000000 + 1 * 10 + 0), MPI_COMM_WORLD, &requests[1]);
                    }
                    if(rank == size - 1)
                    {
                        // send the second row to process size - 2 (the first one is the ghost one, process size - 2 has it already)
                        MPI_Isend(local_solution.data() + N, N, MPI_DOUBLE,
                                size - 2, (1000000 + (size - 1) * 10 + size - 2), MPI_COMM_WORLD, &requests[4*(size - 1) - 2]);

                        // receive the second last row from process size - 2 (that will be placed in the first ghost row of process size - 1)
                        MPI_Irecv(local_solution.data(), N, MPI_DOUBLE,
                                size - 2, (1100000 + (size - 2) * 10 + size - 1), MPI_COMM_WORLD, &requests[4*(size - 1) - 1]);
                    }
                    if(rank != 0 && rank != size - 1)
                    {
                        // send the second last row to process rank + 1 (the last one is the ghost one, process rank + 1 has it already)
                        MPI_Isend(local_solution.data() + local_solution.size() - 2 * N, N, MPI_DOUBLE,
                                rank + 1, (1100000 + rank * 10 + rank + 1), MPI_COMM_WORLD, &requests[4*rank - 2]);
                        // receive the first row from process rank + 1 (that will be placed in the last ghost row of process rank)
                        MPI_Irecv(local_solution.data() + local_solution.size() - N, N, MPI_DOUBLE,
                                rank + 1, (1000000 + (rank + 1) * 10 + rank), MPI_COMM_WORLD, &requests[4*rank -1]);

                        // send the second row to process rank - 1 (the first one is the ghost one, process rank - 1 has it already)
                        MPI_Isend(local_solution.data() + N, N, MPI_DOUBLE,
                                rank - 1, (1000000 + rank * 10 + rank - 1), MPI_COMM_WORLD, &requests[4*rank]);
                        // receive the second last row from process rank - 1 (that will be placed in the first ghost row of process rank)
                        MPI_Irecv(local_solution.data(), N, MPI_DOUBLE,
                                rank - 1, (1100000 + (rank - 1) * 10 + rank), MPI_COMM_WORLD, &requests[4*rank + 1]);
                    }
                    

                    // Wait for all communications to finish.
                    MPI_Waitall(4*(size -1), requests.data(), statuses.data());
                }

                ++n_iterations;
            }

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