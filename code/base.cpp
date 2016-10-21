// TODO
// - times are wrong, too slow
// - temperature printing alignment


// TODO -- Ashish 
// Need to see if initial path needs to be the order in which the nodes were entered OR use the k-nearest neighbour method ( needs extra implementation and more runtime too maybe? )
// 



// tsp.cpp
// Traveling salesman problem with simulated annealing

// usage:  g++ ./tsp.cpp -o tsp
//         ./tsp < input


// input file format:
// a line containing the intger n
// n lines of %d %d, the co-ordinates of each city

// example
// 3        // 3 cities
// 1 1      // co-ords of city 0
// 1 2      // co-ords of city 1
// 3 2      // co-ords of city 2

#include <cstring>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <utility>
#include <algorithm>
#include "omp.h"
//#include "mpi.h"

typedef std::pair<int, int> int_pair;
using namespace std;

// annealing temperature constants
#define TEMP_START    100
#define TEMP_END      0.00001
#define TEMP_COOLING  0.9999

// set to true to use hill climbing instead of annealing
#define HILL_CLIMBING false
#define UNIFORM_V 0.1

// what to print to stdout
// 0: only print the min path cost found
// 1: print inital path cost, min path cost and total time taken
// 2: print all of 1, print each new min path cost and when it is found
// 3: print all of 2, print the current temperature
// NOTE: io is expensive, VERBOSITY 3 considerably slows run time
#define VERBOSITY 3

struct timeval start, end; // track start and end time of program run

class TSP {
    private:
        int      n;                 // num cities
        clock_t  start_time;        // start time of run
        double   temperature;       // starting temperature
        double **distances;         // distance array between cities

        int* cur_path;              // the current path explored
        int* new_path;              // the next path to explore
        int* min_path;              // shortest path found

        void swap(int, int, int[]); // swap cities on a path
        double cost(int[]);        // find the path cost of a path

    public:
        TSP();                      // constructor
       ~TSP();                      // destructor

        double metropolis(float,float,double);          // Metropolis Algorithm -- The Basis of both SA & PT
        //void pa_openMP();
        //void pa_MPI();                  // Find the min path via Parallel Tempering -- According to the method specified in the paper
        void simul_anneal(int,int,int);        // Find the min path via Simulated Annealing -- According to the method specified in the paper
        void anneal();              // find the min path via annealing
        void read_in();             // read in the problem from file

        double distance(int_pair, int_pair);

};

double timedifference_msec(struct timeval start, struct timeval end) {
                // seconds             +  // milliseconds
    return (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
}

// get the path cost of a tsp path
double TSP::cost(int path[]){
    double cost = 0;
    int from_city, to_city;

    for(int i = 0; i < n+1; ++i){
        from_city = path[i];
        to_city   = path[i+1];
        cost += distances[from_city][to_city];
    }
    return cost;
}

// find the shortest path in the tsp (typically using simulated annealing)
void TSP::anneal(){

    // instantiate paths
    cur_path = new int[n+1];
    new_path = new int[n+1];
    min_path = new int[n+1];

    // always use city 0 as the start and end point
    // this does not change throughout the program
    cur_path[0] = cur_path[n] = 0;
    new_path[0] = new_path[n] = 0;
    min_path[0] = min_path[n] = 0;

    // use the order of input as the start implementation
    for(int i = 1; i < n; ++i){
        cur_path[i] = i;
        new_path[i] = i;
        min_path[i] = i;
    }

    // calculate the initial distance of each path
    // (all paths are identical to start)
    float init_cost = cost(cur_path);

    float cur_cost = init_cost;
    float min_cost = init_cost;
    float new_cost = init_cost;

    temperature = TEMP_START;

    // repeatedly try new paths, a higher temperature means that a
    // suboptimal path may be explored in order to escape local maximums
    while(temperature > TEMP_END){

        // pick two random cities to swap (excluding starting / end city)
        int city1 = (rand() % (n-1)) + 1;
        int city2 = (rand() % (n-1)) + 1;

        swap(city1, city2, new_path);

        // calculate new path cost, use new path or don't
        // based on relative costs and temperatures
        new_cost = cost(new_path);
        int gain = new_cost - cur_cost;

        double update_prob = 1 / (1 + pow(M_E, gain/temperature));

        // if new_path is shorter than min the update min
        if(new_cost < min_cost){
            for(int i=1; i < n; ++i){
                min_path[i] = new_path[i];
            }
            min_cost = new_cost;

            if(VERBOSITY >= 2){
                gettimeofday(&end, 0);
                double elapsed = timedifference_msec(start, end);
                printf("\rtime: %.3fs\tpath: %.3f\n", elapsed, min_cost);
            }
        }

        // use to swap between simulated annealing and hill climbing
        bool use_new = HILL_CLIMBING ? new_cost < cur_cost : update_prob > (double) (rand()/(double) RAND_MAX);

        if(use_new){
            // use new_path
            swap(city1, city2, cur_path);
            cur_cost = new_cost;
        } else {
            // revert to cur path
            swap(city1, city2, new_path);
            new_cost = cur_cost;
        }

        if(VERBOSITY >= 3){
            // TODO fix alignment
            int x = 8;
            printf("\r\t\t\t\t\ttemp = %16f", temperature);
        }
        temperature *= TEMP_COOLING;
    }

    if(VERBOSITY >= 1){
        printf("\n\nshortest path: %.3f\n", min_cost);
        printf("initial  path: %.3f\n", init_cost);

        // for(int i = 0; i < n+1; ++i){
        //     printf("%d ", min_path[i]);
        // }

        gettimeofday(&end, 0);
        double elapsed = timedifference_msec(start, end);
        printf("\ntime: %.3fs\n", elapsed);
    } else {
        printf("%.3f\n", min_cost);
    }

}

// swap two cities (a and b) on a path
void TSP::swap(int a, int b, int path[]){
    int temp = path[a];
    path[a] = path[b];
    path[b] = temp;
}

// constructor
TSP::TSP(){
    // because we are comparing runtimes we don't want stochastic behaviour.
    // set static seed
    srand(1001);

    gettimeofday(&start, 0);
}

// destructor
TSP::~TSP(){

}

// read in tsp problem from file
// format described at top of file
void TSP::read_in(){
    scanf("%d", &n);

    // store cities with their coords
    int_pair *cities = new int_pair[n];

    int x,y;

    // read in cities with co-ords
    for(int i = 0; i < n; ++i){
        scanf("%d %d", &x, &y);
        cities[i] = {x, y};
    }

    // instantiate distances array and populate
    distances = new double*[n];
    for(int i=0; i < n; ++i){
        distances[i] = new double[n];
    }

    double d;

    for(int i = 0; i < n; ++i) {
        for(int j = i; j < n; ++j) {
            d = distance(cities[i], cities[j]);

            // graph is undirected, both directions are identical
            distances[i][j] = d;
            distances[j][i] = d;
        }
    }
}

double TSP::distance(int_pair a, int_pair b){
    int x_dist = abs(a.first  - b.first);
    int y_dist = abs(a.second - b.second);

    return sqrt(x_dist * x_dist + y_dist * y_dist);
}


double TSP::metropolis(float new_cost, float old_cost, double temp) {
     // Acceptance probability from old path to new path
    return exp((old_cost - new_cost)/temperature);   // change to min(1,exp((old_cost - new_cost)/temperature));
}

void TSP::simul_anneal(int total_steps, int step_length, int cooling_rate){

    // instantiate paths
    cur_path = new int[n+1];
    new_path = new int[n+1];
    //min_path = new int[n+1];

    // always use city 0 as the start and end point
    // this does not change throughout the program
    cur_path[0] = cur_path[n] = 0;
    new_path[0] = new_path[n] = 0;
    //min_path[0] = min_path[n] = 0;

    // use the order of input as the start implementation
    for(int i = 1; i < n; ++i){
        cur_path[i] = i;
        new_path[i] = i;
       // min_path[i] = i;
    }

    float init_cost = cost(cur_path);
    int step_temp; 
    float cur_cost = init_cost;
    float min_cost = init_cost;
    float new_cost = init_cost;

    double init_temp = 1.5/sqrt(n);  // DEFINE initial temperature according to number of cities
    double curr_temp; 

    step_temp = total_steps/step_length - 1; 

    for (int i=0; i<step_temp ;i++){
        curr_temp = init_temp * pow (cooling_rate, i); 

        for (int j=0; j<=step_length; j++){

            int city1 = (rand() % (n-1)) + 1;
            int city2 = (rand() % (n-1)) + 1;
            swap(city1, city2, cur_path);

            new_cost = cost(cur_path);
            
            if(new_cost < min_cost){
                for(int i=1; i < n; ++i){
                    min_path[i] = cur_path[i];
                }
                min_cost = new_cost;
            }

            if (metropolis(new_cost,cur_cost,curr_temp) > UNIFORM_V){
                //cur_path = new_path;
                cur_cost = new_cost; 
            }
            else swap(city1,city2,cur_path);
        }

    }

    // Didn't go into implementing the VERBOSITY prints needed, stuck to  -- VERBOSITY = 1

    printf("\n\nshortest path: %.3f\n", min_cost);
    printf("initial  path: %.3f\n", init_cost);


}


void TSP::pa_openMP(int total_steps, int swap_length, int cooling_rate){

    cur_path = new int[n+1];
    new_path = new int[n+1];
    //min_path = new int[n+1];

    // always use city 0 as the start and end point
    // this does not change throughout the program
    cur_path[0] = cur_path[n] = 0;
    new_path[0] = new_path[n] = 0;
        //min_path[0] = min_path[n] = 0;

    // use the order of input as the start implementation
    for(int i = 1; i < n; ++i){
        cur_path[i] = i;
        new_path[i] = i;
        //min_path[i] = i;
    }

    float init_cost = cost(cur_path);
    int step_temp; 
    float cur_cost = init_cost;
    float min_cost = init_cost;
    float new_cost = init_cost;

    double init_temp = 1.5/sqrt(n);  // DEFINE initial temperature according to number of cities
    double curr_temp; 

    step_temp = total_steps/swap_length - 1; 
    int temp_array[5] = {0.0025, 0.004, 0.006, 0.008,0.01}; 

    for (int i=0; i<step_temp ;i++){
        //vector<vector<int>> distance_array(n, vector<int>(n+1));
        curr_temp = init_temp * pow (cooling_rate, i); 
        //#pragma omp parallel for firstprivate (curr_temp,cur_path,new_path) lastprivate (cur_path) (5)
        for (j=0; j<=swap_length; j++){

            int city1 = (rand() % (n-1)) + 1;
            int city2 = (rand() % (n-1)) + 1;
            swap(city1, city2, cur_path);

            new_cost = cost(new_path);

            if (metropolis(new_cost,cur_cost,curr_temp) > UNIFORM_V){
            //    cur_path = new_path;
                cur_cost = new_cost; 
            }
            else swap(city1,city2,cur_path);
            // HERE,  AFTER the j-loop finishes execution, we need to exchange the path of any two RANDOM paths
        }
            // TO DO TO DO 
    }

}*/

int main(){
    TSP tsp;
    tsp.read_in();
    tsp.anneal();
    tsp.simul_anneal(1000000,160000,0.95);
    return 0;
}
