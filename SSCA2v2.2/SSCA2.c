#include "defs.h"

#ifdef _OPENMP
int NUM_THREADS;
#endif

int main(int argc, char** argv)
{
    /* Data structure for storing generated tuples in the 
     * Scalable Data Generation Stage -- see defs.h */
    graphSDG* SDGdata;  

    /* The graph data structure -- see defs.h */
    graph* G;
  
    /* Kernel 2 output */
    edge *maxIntWtList;
    INT_T maxIntWtListSize;

    /* Kernel 4 output */
    DOUBLE_T *BC;
    
    DOUBLE_T elapsed_time;

#ifdef _OPENMP
    if (argc != 3) {
        
        fprintf(stderr, "Usage: ./SSCA2 <No. of threads> <SCALE>\n");
        exit(-1);
    }
    NUM_THREADS = atoi(argv[1]);
    SCALE = atoi(argv[2]);
    omp_set_num_threads(NUM_THREADS);
#else
    if (argc != 2) {
        fprintf(stderr, "Usage: ./SSCA2 <SCALE>\n");
        exit(-1);
    }
    SCALE = atoi(argv[1]);
#endif

    /* ------------------------------------ */
    /*  Initialization -- Untimed           */
    /* ------------------------------------ */

    fprintf(stderr, "\nHPCS SSCA Graph Analysis Benchmark v2.2\n");
    fprintf(stderr, "Running...\n\n");

    init(SCALE);

#ifdef _OPENMP
    fprintf(stderr, "\nNo. of threads: %d\n", NUM_THREADS);
#endif
    fprintf(stderr, "SCALE: %d\n\n", SCALE);
 
    /* -------------------------------------------- */
    /*  Scalable Data Generator -- Untimed          */
    /* -------------------------------------------- */

#ifndef VERIFYK4
    fprintf(stderr, "Scalable Data Generator -- ");
    fprintf(stderr, "genScalData() beginning execution...\n");

    SDGdata  = (graphSDG *) malloc(sizeof(graphSDG));
    elapsed_time = genScalData(SDGdata);

    fprintf(stderr, "\n\tgenScalData() completed execution\n");
    fprintf(stderr, 
            "\nTime taken for Scalable Data Generation is %9.6lf sec.\n\n", 
            elapsed_time);
#else
    fprintf(stderr, "Generating 2D torus for Kernel 4 validation -- ");
    fprintf(stderr, "gen2DTorus() beginning execution...\n");

    SDGdata = (graphSDG *) malloc(sizeof(graphSDG));
    elapsed_time = gen2DTorus(SDGdata);

    fprintf(stderr, "\n\tgen2DTorus() completed execution\n");
    fprintf(stderr, 
            "\nTime taken for 2D torus generation is %9.6lf sec.\n\n", 
            elapsed_time);
#endif
    
    
    /* ------------------------------------ */
    /*  Kernel 1 - Graph Construction       */
    /* ------------------------------------ */
  
    /* From the SDG data, construct the graph 'G'  */
    fprintf(stderr, "\nKernel 1 -- computeGraph() beginning execution...\n");

    G = (graph *) malloc(sizeof(graph));
    /* Store the SDG edge lists in a compact representation 
     * which isn't modified in subsequent Kernels */
    elapsed_time = computeGraph(G, SDGdata);

    fprintf(stderr, "\n\tcomputeGraph() completed execution\n");
    fprintf(stderr, "\nTime taken for Kernel 1 is %9.6lf sec.\n\n", 
            elapsed_time);
    
    free(SDGdata);

    /* ---------------------------------------------------- */
    /*  Kernel 2 - Find max edge weight                     */
    /* ---------------------------------------------------- */
  
    fprintf(stderr, "\nKernel 2 -- getStartLists() beginning execution...\n");
  
    /* Initialize vars and allocate temp. memory for the edge list */
    maxIntWtListSize = 0;
    maxIntWtList = (edge *) malloc(sizeof(edge));

    elapsed_time = getStartLists(G, &maxIntWtList, &maxIntWtListSize);

    fprintf(stderr, "\n\tgetStartLists() completed execution\n\n");
    fprintf(stderr, "Max. int wt. list size is %d\n", maxIntWtListSize);
    fprintf(stderr, "\nTime taken for Kernel 2 is %9.6lf sec.\n\n", 
            elapsed_time);

    /* ------------------------------------ */
    /*  Kernel 3 - Graph Extraction         */
    /* ------------------------------------ */
  
    fprintf(stderr, "\nKernel 3 -- findSubGraphs() beginning execution...\n");

    elapsed_time = findSubGraphs(G, maxIntWtList, maxIntWtListSize);

    fprintf(stderr, "\n\tfindSubGraphs() completed execution\n");
    fprintf(stderr, "\nTime taken for Kernel 3 is %9.6lf sec.\n\n", 
            elapsed_time);
     
    free(maxIntWtList);
       
    /* ------------------------------------------ */
    /*  Kernel 4 - Betweenness Centrality         */
    /* ------------------------------------------ */
  
    fprintf(stderr, "\nKernel 4 -- betweennessCentrality() "
            "beginning execution...\n");  
    
    BC = (DOUBLE_T *) calloc(N, sizeof(DOUBLE_T));
    elapsed_time = betweennessCentrality(G, BC);
    fprintf(stderr, "\nTime taken for Kernel 4 is %9.6f sec.\n\n", 
            elapsed_time);
    fprintf(stderr, "TEPS score for Kernel 4 is %lf\n\n", 
            7*N*((long)1<<K4approx)/elapsed_time);

    /***************************************************************
     * Validation LN and PHJK (execution times on 4-threads i7-4790)
     * This is a table of reference outputs, but since the algorithm
     * is non-deterministic they are only representative.
     * You will need to add to this table if you run for larger
     * SCALE.
     **************************************************************/
    int i = 0;
    double maxBC = BC[0], minBC = BC[0], sum = 0., avg = 0.;
    long int intVal = 0., intMaxBC = 0, intMinBC = 0, intAvg = 0;
    #define NResults 50
    long int prestoredResults[NResults][3] = { 
    /* SCALE=1 */     { 0, 0, 0 }, 
    /* SCALE=2 */     { 0, 2500, 1250 }, 
    /* SCALE=3 */     { 0, 30999, 6249 }, 
    /* SCALE=4 */     { 0, 48028, 14812 }, 
    /* SCALE=5 */     { 0, 133238, 44000 }, 
    /* SCALE=6 */     { 0, 368702, 115578 }, 
    /* SCALE=7 */     { 0, 1790304, 257476 }, 
    /* SCALE=8 */     { 0, 5326002, 547863 }, 
    /* SCALE=9 */     { 0, 21598987, 1237150 }, 
    /* SCALE=10 */     { 0, 40382717, 2415442 }, 
    /* SCALE=11 */     { 0, 77220870, 2828273 }, 
    /* SCALE=12 */     { 0, 155192044, 3002415}, 
    /* SCALE=13 */     { 0, 237379503, 3098723}, 
    /* SCALE=14 */     { 0, 375125133, 3250859 }, 
    /* SCALE=15 */     { 0, 701218499, 3277281 }, 
    /* SCALE=16 */     { 0, 1158382011, 3338237 }, 
    /* SCALE=17 */     { 0, 1791767763, 3408419 }, 
    /* SCALE=18 */     { 0, 3341777600, 3491043 }, 
    /* SCALE=19 */     { 0, 5511190429, 3522221  }, /* ca.100 seconds */
    /* SCALE=19 */     /* { 0, 5321333625, 3523534  }, another run */
    /* SCALE=20 */     { 0, 11269213138, 3577907 }, /* ca.200 seconds */
    /* SCALE=20 */     /* { 0, 10994335868, 3553423 }, another run */
    /* SCALE=21 */     { 0, 15593566622, 3556022 } /* 2000 seconds */
    };
    
    fprintf(stderr, "N=%d, SCALE=%d\n", N, SCALE); 
    for(i=0; i<N; i++)
    {
	if(BC[i]>maxBC)
	{
	    maxBC = BC[i];
	}
	if(BC[i]<minBC)
	{
	    minBC = BC[i];
	}
	sum += BC[i];
    }    
    avg = sum / N; 

    // It is easier to convert them in integer to do the comparison and store the array in prestoredResults.
    intMaxBC = maxBC * 1000;
    intMinBC = minBC * 1000;
    intAvg = avg * 1000;

    fprintf(stderr, "min=%ld, max=%ld, avg=%ld\n", intMinBC, intMaxBC, intAvg); 
    fprintf(stderr, "This result should match:\n"); 
    fprintf(stderr, "min=%ld, max=%ld, avg=%ld\n\n", 
	prestoredResults[SCALE-1][0], prestoredResults[SCALE-1][1], prestoredResults[SCALE-1][2]); 
    fprintf(stderr, "%ld, %ld, %ld\n", intMinBC, intMaxBC, intAvg); 

    /* PHJK: the betweenness centrality algorithm is slightly non-deterministic
     * so we check that the min, max and avg are within reasonable bounds.
     * This might need adjusting if you manage very large problems.
     */
    if(labs(intMinBC - prestoredResults[SCALE-1][0]) <= 1+(intMinBC/10) && 
       labs(intMaxBC - prestoredResults[SCALE-1][1]) <= 1+(intMaxBC/5) && 
       labs(intAvg - prestoredResults[SCALE-1][2]) <= 1+(intAvg/20))
    {
        fprintf(stderr, "Kernel 4 validation successful!\n");
    } else {
        fprintf(stderr, "Kernel 4 failed validation!\n");
    }
 
    
    free(BC);
    
    /* -------------------------------------------------------------------- */
    
    /* free(G->numEdges);  PHJK: crashes under simplescalar
    free(G->endV); 
    free(G->weight);
    free(G); */

    return 0;
}
