#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include "util.h"
#include <stdbool.h>

// mpi header file
#include <mpi.h>

// sssp_parallel.c program an example call is
// mpirun -np 4 ./sssp_parallel ANY_INPUT_FILENAME.txt 4 ANY_OUTPUT_FILENAME.txt
// yes the example on assignment does not call via mpirun, but it is wrong

// main function
int main(int argc, char *argv[])
{

    // check if the number of arguments is correct
    if (argc != 4)
    {
        printf("Usage: %s <input_file> <source_vertex> <output_file>\n", argv[0]);
        return 1;
    }

    // MPI variables
    int rank, size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes

    int sourceVertex = atoi(argv[2]);

    // start timer
    double startTime = MPI_Wtime();

    // offsets is the rows
    // edges is the columns
    // why did they named it like that? To confuse us?
    int numVertices, numEdges, *offsets, *edges, *weights;

    if (rank == 0)
    {
        printf("sourceVertex = %d\n", sourceVertex);
        // read the input file
        int success = read_file(argv[1], &numVertices, &numEdges, &offsets);
        if (!success)
        {
            printf("Error reading file %s\n", argv[1]);
            return 1;
        }
        else
        {
            edges = &(offsets[numVertices + 1]);
            weights = &(offsets[numVertices + 1 + numEdges]);
            // print the edges, weights and offsets arrays
            // printf("\noffsets: ");
            for (int i = 0; i < numVertices + 1; i++)
            {
                // printf("%3d ", offsets[i]);
            }
            // printf("\nedges:   ");
            for (int i = 0; i < numEdges; i++)
            {
                // printf("%3d ", edges[i]);
            }
            // printf("\nweights: ");
            for (int i = 0; i < numEdges; i++)
            {
                // printf("%3d ", weights[i]);
            }
            // printf("\n\n");
        }

        // print each node and its neighbors
        // for (int i = 0; i < numVertices; i++)
        // {
        //     printf("Node %d: ", i);
        //     for (int j = offsets[i]; j < offsets[i + 1]; j++)
        //     {
        //         printf("(%d, %d)", edges[j], weights[j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n\n");
    } // end of rank 0 reading the file

    // Broadcast the number of vertices, edges and offsets to all processes
    MPI_Bcast(&numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // we assume that the number of vertices is divisible by the number of processes
    int numVerticesPerProcess = numVertices / size;

    // printf("rank %3d: numVerticesPerProcess = %d/%d=%3d numedges=%3d\n", rank, numVertices, size, numVerticesPerProcess, numEdges);
    // MPI_Barrier(MPI_COMM_WORLD);

    // now create displs and sendcounts for scatterv
    int *sendcounts = (int *)malloc(size * sizeof(int));
    int *displs = (int *)malloc(size * sizeof(int));
    int *offsetsendcounts = (int *)malloc(size * sizeof(int));
    int *offsetsdispls = (int *)malloc(size * sizeof(int));
    if (rank == 0)
    {
        // initialize all to 0
        // printf("%d: initialize all to 0\n", rank);
        for (int i = 0; i < size; i++)
        {
            sendcounts[i] = 0;
            displs[i] = 0;
            offsetsendcounts[i] = 0;
            offsetsdispls[i] = 0;
        }

        // using the offsets array, we can calculate the sendcounts and displs
        // printf("%d: compute sendcounts etc\n", rank);
        for (int i = 0; i < size; i++)
        {
            // printf("%d: offset[%d]=%3d\n", rank, i * numVerticesPerProcess, offsets[i * numVerticesPerProcess]);

            sendcounts[i] = offsets[(i + 1) * numVerticesPerProcess] - offsets[i * numVerticesPerProcess];
            displs[i] = offsets[i * numVerticesPerProcess];
            // printf("rank %3d:       sendcounts[%3d] = %3d        displs[%3d] = %3d \n", rank, i, sendcounts[i], i, displs[i]);

            // +1 is due "Remember that rows array has an additional index that array points
            // to one plus last index of columns array."
            offsetsendcounts[i] = numVerticesPerProcess + 1;
            offsetsdispls[i] = i * numVerticesPerProcess;
            // printf("rank %3d: offsetsendcounts[%3d] = %3d offsetsdispls[%3d] = %3d \n\n", rank, i, offsetsendcounts[i], i, offsetsdispls[i]);
        }
    } // all sendcounts etc. are calculated by rank 0

    // now broadcast the sendcounts and displs to all processes
    MPI_Bcast(sendcounts, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offsetsendcounts, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offsetsdispls, size, MPI_INT, 0, MPI_COMM_WORLD);

    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("rank %3d:       sendcounts[%3d] = %3d displs[%3d] = %3d \n", rank, rank, sendcounts[rank], rank, displs[rank]);
    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("rank %3d: offsetsendcounts[%3d] = %3d offsetsdispls[%3d] = %3d \n", rank, rank, offsetsendcounts[rank], rank, offsetsdispls[rank]);
    // MPI_Barrier(MPI_COMM_WORLD);

    // now scatter the offsets, edges and weights arrays

    int *offsetsPerProcess = (int *)malloc((offsetsendcounts[rank]) * sizeof(int));
    int *edgesPerProcess = (int *)malloc(sendcounts[rank] * sizeof(int));
    int *weightsPerProcess = (int *)malloc(sendcounts[rank] * sizeof(int));

    MPI_Scatterv(offsets, offsetsendcounts, offsetsdispls, MPI_INT, offsetsPerProcess, offsetsendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(edges, sendcounts, displs, MPI_INT, edgesPerProcess, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(weights, sendcounts, displs, MPI_INT, weightsPerProcess, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

    // by now each process has its own offsetsPerProcess, edgesPerProcess and weightsPerProcess
    // to help me understand, I will print them
    // todo comment out
    // printf("rank %3d has %d offsetsPerProcess: \n", rank, offsetsendcounts[rank]);
    // MPI_Barrier(MPI_COMM_WORLD); // all should wait TODO debug

    // for (int i = 0; i < offsetsendcounts[rank]; i++)
    // {
    //     printf("%5d:%3d offset %3d \n", rank, i, offsetsPerProcess[i]);
    // }

    // MPI_Barrier(MPI_COMM_WORLD); // all should wait TODO debug

    // for (int i = 0; i < sendcounts[rank]; i++)
    // {
    //     printf("%3d: i=%3d to %3d w=%3d \n", rank, i + offsetsPerProcess[0], edgesPerProcess[i], weightsPerProcess[i]);
    // }
    // printf("\n");

    // At  this point each process has its own offsetsPerProcess, edgesPerProcess and weightsPerProcess
    // now we will use similar code as in sequential to calculate the shortest path
    // The description of the algorithm is in the assignment:
    //
    // "Each processors will be responsible for iterative computation of the shortest distances in the set of vertices they
    // have. In doing so, a global distance vector D of size |V | should be initialized in each processor as
    // in the serial algorithm where only D[sourceVertex] equals to 0. For each iteration, each processor
    // computes a local distance vector R of size |V |/p where p is the number of processors. R gives the
    // resulting distance values of the a processor’s vertices at the current iteration. Before moving on to
    // next iteration, we should use all-gather operation to update the global distance vector D. In this
    // operation each processor sends their local R vectors to all other processors and receives the local
    // R vectors of all other processors with placing them onto the global distance vector D."
    int *D = (int *)malloc(numVertices * sizeof(int));
    int *R = (int *)malloc(numVerticesPerProcess * sizeof(int));
    // initialize D
    for (int i = 0; i < numVertices; i++)
    {
        D[i] = INT_MAX;
    }
    D[sourceVertex] = 0;

    if (rank == 0)
    {
        // printf("D: ");
        for (int i = 0; i < numVertices; i++)
        {
            // printf("D[%d]=%3d \n", i, D[i]);
        }
        // printf("\n");
    }
    // initialize R
    for (int i = 0; i < numVerticesPerProcess; i++)
    {
        R[i] = INT_MAX;
        if (i + rank * numVerticesPerProcess == sourceVertex)
        {
            R[i] = 0; // to make sure that the source vertex is 0
        }
    }

    // now we will start the iterations
    int *changed = malloc(sizeof(int));
    *changed = 1;
    for (int i = 0; i < numVertices - 1 && *changed > 0; i++)
    {
        *changed = 0;

        // // update R from D
        // for (int j = 0; j < numVerticesPerProcess; j++)
        // {
        //     R[j] = D[j + offsetsPerProcess[0]];
        // }

        // for each vertex in the process
        for (int j = 0; j < numVerticesPerProcess; j++)
        {
            int l = offsetsPerProcess[j];
            int r = offsetsPerProcess[j + 1];

            // for each edge of the vertex
            for (int k = l; k < r; k++)
            {
                int u = edgesPerProcess[k - offsetsPerProcess[0]];
                int w = weightsPerProcess[k - offsetsPerProcess[0]];

                // Either j or u supposed to be read from D but not sure which rn
                if (D[u] + w > 0 && D[u] + w < R[j])
                {
                    // log the change
                    *changed = 1;
                    // printf("%d changed distance of %d from %d to %d by coming from %d. changed=%d \n", rank, j + numVerticesPerProcess * rank, R[j], D[u] + w, u, *changed);
                    R[j] = D[u] + w;
                }
                // // if (R[j] + w >= 0 && R[j] + w < R[u - offsetsPerProcess[0]])
                // if (R[j] + w >= 0 && R[j] + w < D[u])
                // {
                //     // log the change
                //     printf("%d changed distance of %d from %d to %d \n", rank, u, R[u - offsetsPerProcess[0]], R[j] + w);
                //     // printf("%d changed distance of %d from %d to %d \n", rank, u, D[u], R[j] + w);
                //     R[u - offsetsPerProcess[0]] = R[j] + w;
                //     // D[u] = R[j] + w;
                //     changed = true;
                // }
            }
        }
        // printf("%d: at the end of iteration %d R: ", rank, i);

        for (int j = 0; j < numVerticesPerProcess; j++)
        {
            // printf("D[%d]=R[%d]=%d, ", j + rank * numVerticesPerProcess, j, R[j]);
        }
        // printf("\n");

        // update global distance D with Allgather
        MPI_Allgather(R, numVerticesPerProcess, MPI_INT, D, numVerticesPerProcess, MPI_INT, MPI_COMM_WORLD);

        // notify all processes that D has been updated
        // that is AllReduce the changed variable
        MPI_Allreduce(MPI_IN_PLACE, changed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        //  printf("%d: at the end of iteration %d changed=%d ", rank, i, *changed);
    }

    // now we will print the results
    if (rank == 0)
    {
        printResults(argv[3], D, numVertices);
    }
    // end timer
    double endTime = MPI_Wtime();
    double totalTime = endTime - startTime;
    printf("rank %3d: totalTime = %f\n", rank, totalTime);
    MPI_Finalize();
} // end main