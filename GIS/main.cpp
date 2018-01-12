#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <fstream>
#include <string>
#include <limits>
using namespace std;


#define NUM_INTERATION 10000	// maximal number of iterations
#define TIME_TRY 200			// maximal number of iterations with no improvement

int N; // number of cities
double **costMatrix; // cost of each edge
int **tabu_list;
double score=0.0, bestSolverScore=0.0;
int *v, *foundSolution; // vertices, found path
double infinity = 1e+38;


string generateData(int n, double percent, int min=0, int max=100){
	string filename = "dane" + std::to_string(n) + "_" + std::to_string((int)percent)+".txt";


	double **cost = new double*[n];	// n - number of cities, create n x n cost matrix
	for(int i = 0; i < n; ++i){
		cost[i] = new double[n];
	}

	ofstream myfile;
	myfile.open(filename.c_str());
	myfile << n << "\n";

	srand( time( NULL ) );
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			cost[i][j] = (max - min) * ( (double)rand() / (double)RAND_MAX ) + min; // random cost
		}
	}

	if(percent>1) percent = percent/100; // % of one-way edges in G
	int tabuNum = percent*n*(n-1) + 1; // number of one-way edges
	int x=0, y=0, z=0;
	while(z!=tabuNum){			// create one-way edges randomly (non overlapping)
		x = (n-1) * ( (double)rand() / (double)RAND_MAX );
		y = (n-1) * ( (double)rand() / (double)RAND_MAX );
		if(cost[x][y]!=infinity && x!=y && y!=(x+1)){
			cost[x][y] = infinity;
			++z;
		}
	}

	for(int i=0; i<n; ++i){ // cycle required
		cost[i][i] = 0;
		cost[i][i+1] = 1;
	}
	cost[n-1][0]=1;

	for(int i=0; i<n; ++i){ // write to file
		for(int j=0; j<n; ++j){
			myfile <<cost[i][j] << " ";
		}
		myfile<<"\n";
	}
	myfile.close();

	for(int i = 0; i < n; ++i){
			delete[] cost[i];
	}
	delete[] cost;

	return filename;
}


void readCostMatrix(string fileName){
    ifstream file(fileName.c_str());

    if(!file.is_open()){
        perror ("The following error occurred");
        exit(0);
    }

    if(file.is_open()){
    	file >> N; // known number of cities

    	costMatrix = new double*[N];
    	for(int i = 0; i < N; ++i){		// dynamically allocate cost matrix
    		costMatrix[i] = new double[N];
    	}

    	v = new int[N];
    	foundSolution = new int[N];
    	score = 0;

        for(int i = 0; i < N; ++i){		// read file
        	for(int j = 0; j<N; ++j){
        		file >> costMatrix[i][j];
        	}
        }
    }

}

double getScore(int *v){ // cost of given cycle
	double score = 0;
	for(int i = 0; i < (N - 1); ++i){		// sum cost of each edge in cycle
		score += costMatrix[v[i]][v[i+1]];
	}
	return score += costMatrix[v[N-1]][v[0]];
}

void resetTabuList(){ // empty tabu list
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < N; ++j){
			tabu_list[i][j] = 0;
		}
	}
}

void initTabuList(){ // dynamically allocate tabu list
	tabu_list = new int*[N];

	for(int i = 0; i < N; ++i){
		tabu_list[i] = new int[N];
	}
	resetTabuList();
}

void initSolution(){ // initial (random) path
	for(int i = 0; i < N; i++){
		v[i] = i;
	}
	srand(time(NULL));
	for(int i = (N - 1); i >= 0; --i){
		int j = rand() % N;
		swap(v[i], v[j]);
	}

	score = getScore(v);
}


int* getBestNearbySolution(int it){ // search for neighbours
	int *v_temp = new int[N]; // copy of current solution

	for(int i = 0; i < N; ++i){
		v_temp[i]=v[i];
	}

	double bestScore = std::numeric_limits<double>::max(); // worst case scenario
	int vertexA = 0, vertexB = 1;	// vertices to swap

	for(int i = 0; i < N; ++i){
		for(int j = (i+1); j < N; ++j){ // directed graph, swap of adjacent edges possible

			swap(v_temp[i], v_temp[j]); //swap for new solution

			double currentScore = getScore(v_temp);

			// found solution is better and not tabu or the best of all found
			if( (bestScore > currentScore && tabu_list[i][j] <= it) || currentScore < bestSolverScore){
				vertexA = i;	// remember best neighbour
				vertexB = j;
				bestScore = currentScore;
			}
			swap(v_temp[j], v_temp[i]); // back to original solution
		}
	}

	tabu_list[vertexA][vertexB] = (it + 3*N); 	// update tabu list - 3N = tabu length
	swap(v_temp[vertexA], v_temp[vertexB]);		// swap for best neighbour

	return v_temp;
}

double solveTSP(int numCandidate){

	int *v_temp = new int[N];
	double bestSolutionScore = getScore(v);

	for(int loopCount = 0; loopCount < numCandidate; ++loopCount){
		initSolution();
		resetTabuList();

		int countTime = 0; // count times that solver solution is not improved
		bestSolverScore = std::numeric_limits<double>::max();

		for(int i = 0; i < NUM_INTERATION; ++i){

			v_temp = getBestNearbySolution(i);

			double score = getScore(v_temp);

			if(score < bestSolverScore){ // local solution better
				bestSolverScore = score;
				countTime = 0;

				for(int j = 0; j < N; ++j){		// update local solution
					v[j]=v_temp[j];
				}

				if(bestSolverScore < bestSolutionScore){ // found solution is better than  global solution
					for(int j = 0; j < N; ++j){
						foundSolution[j] = v[j];	// update global solution
					}
					bestSolutionScore = bestSolverScore;	// update objective function
				}
			}else{ // no improvement
				++countTime;
				if(countTime > TIME_TRY){
					cout<<"countTime: "<<i<<endl;
					break;
				}
			}
		}

	}
	return bestSolutionScore;
}



int main(int argc, char* argv[]){
	// tabu search for generated costMatrix
	string fn = generateData(100, 0);
	readCostMatrix(fn);

	//tabu search for existing file
	//readCostMatrix("dane50_0.txt");


	initSolution();
	cout<<"Pocz¹tkowe rozwi¹zanie: ";
	for(int i = 0; i < N; ++i){
		cout<<v[i]<<" ";
	}
	cout<<endl<<"Koszt trasy pocz¹kowej: "<<getScore(v)<<endl;

	initTabuList();

	time_t start, end;
	const clock_t begin_time = clock();
	time(&start);

	double best = solveTSP(5);		// number of tries

	time(&end);
	cout << "Czas wykonania: " << difftime(end,start)<<endl;
	cout << "Czas procesora: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC <<endl;

	cout<<"Koñcowe rozwi¹zanie: ";
	for(int i = 0; i < N; ++i){
		cout<<foundSolution[i]<<" ";
	}
	cout<<endl<<"Koszt trasy: "<<best<<endl;

	return 0;
}
