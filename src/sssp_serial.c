#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "util.h"	
#include <stdbool.h>

    
// an example call is
// # ./sssp_serial ANY_INPUT_FILENAME.txt 4 ANY_OUTPUT_FILENAME.txt

// main function
int main(int argc, char *argv[]) {
    // check if the number of arguments is correct
    if (argc != 4) {
        printf("Usage: %s <input_file> <source_vertex> <output_file>\n", argv[0]);
        return 1;
    }
    
    int sourceVertex = atoi(argv[2]);
    
    // offsets is the rows
    // edges is the columns
    // why did they named it like that? To confuse us?

    // read the input file
    int numVertices, numEdges, *offsets, *edges, *weights;
    int success = read_file(argv[1], &numVertices, &numEdges, &offsets);
    if (!success) {
        printf("Error reading file %s\n", argv[1]);
        return 1;
    } else {
        edges = &(offsets[numVertices + 1]);
        weights = &(offsets[numVertices + 1 + numEdges]);
        // print the edges, weights and offsets arrays
        printf("offsets: ");
        for (int i = 0; i < numVertices + 1; i++) {
            printf("%3d ", offsets[i]);
        }
        printf("\nedges:   ");
        for (int i = 0; i < numEdges; i++) {
            printf("%3d ", edges[i]);
        }
        printf("\nweights: ");
        for (int i = 0; i < numEdges; i++) {
            printf("%3d ", weights[i]);
        }
        printf("\n");

    }

    // print each node and its neighbors
    for (int i = 0; i < numVertices; i++) {
        printf("Node %d: ", i);
        for (int j = offsets[i]; j < offsets[i + 1]; j++) {
            printf("(%d, %d)", edges[j], weights[j]);
        }
        printf("\n");
    }

    // allocate memory for the distances array
    int *distances = (int *) malloc(numVertices * sizeof(int));
    if (distances == NULL) {
        printf("Error allocating memory for distances\n");
        return 1;
    }

    // initialize the distances array
    for (int i = 0; i < numVertices; i++) {
        distances[i] = INT_MAX;
    }


    // start the timer
    clock_t start = clock();

    // run the algorithm
    // todo move to another method
    distances[sourceVertex] = 0; // D[sourceVertex] <-0
    bool changed = true;
    
    for(int i=0; i < numVertices - 1 && changed  ;i++){ // repeat |V|−1 times
        changed = false;
        for(int j=0; j < numVertices; j++){ // for each vertex index i in rows do
            int l = offsets[j];
            int r = offsets[j+1];
            for(int k=l; k < r; k++){ // find the minimum distance
                int u = edges[k];
                int w = weights[k];
                if( // distances[j] != INT_MAX 
                    distances[j] + w >= 0
                    && distances[u] > distances[j] + w){
                    // log the change
                    printf("Changed distance of %d from %d to %d due w=%d and d..[%d]=%d\n", u, distances[u], distances[j] + w, w, j, distances[j]);
                    distances[u] = distances[j] + w;
                    changed = true;
                }
            } 
        }

    }

    // stop the timer
    clock_t end = clock();

    // print the results
    printResults(argv[3], distances, numVertices);

    // print the time
    printf("Time: %f ms\n", (double) (end - start) / CLOCKS_PER_SEC * 1000);

    // free the memory
    free(distances);

    return 0;
}
