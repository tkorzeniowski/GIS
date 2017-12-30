#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <fstream>
#include <string>
#include <limits>
using namespace std;


#define TABU_LENGTH 30
#define NUM_INTERATION 10000
//#define PENAL_LONG_TERM 10
#define LONG_TERM_LENGTH 100
#define TIME_TRY 2000

int N; // number of cities
double **costMatrix; // cost of each edge
int **tabu_list;//, **tabu_f_list;
double score=0.0, bestSolverScore=0.0;
int *v, *foundSolution; // vertices, found path
double infinity = 1e+38;


string generateData(int n, double percent, int min=0, int max=100){
	string filename = "dane" + std::to_string(n) + "_" + std::to_string((int)percent)+".txt";


	double **cost = new double*[n];
	for(int i = 0; i < n; ++i){
		cost[i] = new double[n];
	}

	ofstream myfile;
	myfile.open(filename.c_str());
	myfile << n<<"\n";

	srand( time( NULL ) );
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			cost[i][j] = (max - min) * ( (double)rand() / (double)RAND_MAX ) + min; // random cost
		}
	}

	if(percent>1) percent = percent/100; // % of one-way edges in G
	int tabuNum = percent*n*(n-1) + 1; // number of one-way edges
	int x=0, y=0, z=0;
	while(z!=tabuNum){
		x = (n-1) * ( (double)rand() / (double)RAND_MAX );
		y = (n-1) * ( (double)rand() / (double)RAND_MAX );
		if(cost[x][y]!=infinity && x!=y && y!=(x+1)){
			cost[x][y]=infinity;
			++z;
		}
	}

	for(int i=0; i<n; ++i){ // cycle required
		cost[i][i]=0;
		cost[i][i+1]=1;
	}
	cost[n-1][0]=1;

	for(int i=0; i<n; ++i){ // write to file
		for(int j=0; j<n; ++j){
			myfile <<cost[i][j] << " ";
		}
		myfile<<"\n";
	}
	myfile.close();

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
    	for(int i = 0; i < N; ++i){
    		costMatrix[i] = new double[N];
    	}

    	v = new int[N];
    	foundSolution = new int[N];
    	score = 0;

        for(int i = 0; i < N; ++i){
        	for(int j = 0; j<N; ++j){
        		file >> costMatrix[i][j];
        	}
        }
    }

}

double getScore(int *v){ // cost of given cycle
	double score = 0;
	for(int i = 0; i < (N - 1); ++i){
		score += costMatrix[v[i]][v[i+1]];
	}
	return score += costMatrix[v[N-1]][v[0]];
}

void resetTabuList(){ // empty tabu list
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < N; ++j){
			tabu_list[i][j] = 0;
			//tabu_f_list[i][j] = 0;
		}
	}
}

void initTabuList(){ // dynamically allocate tabu list
	tabu_list = new int*[N];
	//tabu_f_list = new int*[N];
	for(int i = 0; i < N; ++i){
		//tabu_f_list[i] = new int[N];
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


//int* getBestNearbySolution(int* v, int it){ // search for neighbours
int* getBestNearbySolution(int it){ // search for neighbours
	int *v_temp = new int[N]; // copy of current solution

	for(int i = 0; i < N; ++i){
		v_temp[i]=v[i];
	}

	double bestScore = std::numeric_limits<double>::max();
	int vertexA = 0, vertexB = 1;

	for(int i = 0; i < N; ++i){
		for(int j = (i+1); j < N; ++j){ // directed graph, swap adjacent edges possible

			swap(v_temp[i], v_temp[j]); //swap for new solution

			double currentScore = getScore(v_temp);
			double penalScore = currentScore ;//+ PENAL_LONG_TERM * tabu_f_list[i][j];
			// found solution is better and not tabu or the best of all
			if( (bestScore > penalScore && tabu_list[i][j] <= it) || currentScore < bestSolverScore){
				vertexA = i;
				vertexB = j;
				bestScore = penalScore;
				//tabu_list[i][j] = (it + TABU_LENGTH);
				//tabu_list[j][i] = (it + TABU_LENGTH);

				tabu_list[i][j] = (it + 3*N);
				//tabu_list[j][i] = (it + 3*N);
			}
			swap(v_temp[j], v_temp[i]); // back to original solution
			//if(tabu_f_list[i][j] > 0 && it > LONG_TERM_LENGTH) tabu_f_list[i][j] -= 1;
		}
	}
	//tabu_f_list[vertexA][vertexB] += 2;
	swap(v_temp[vertexA], v_temp[vertexB]);//s->swapSolve( vertexA, vertexB );
	return v_temp;//s;
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
			//v_temp = getBestNearbySolution(v, i);
			v_temp = getBestNearbySolution(i);

			double score = getScore(v_temp);

			if(score < bestSolverScore){ // local solution better
				bestSolverScore = score;
				countTime = 0;

				if(bestSolverScore < bestSolutionScore){ // found solution is better than  global solution
					for(int j = 0; j < N; ++j){
						v[j]=v_temp[j];//bestSolution.set(j,s->getV(j));
						foundSolution[j] = v[j];
					}
					bestSolutionScore = bestSolverScore;
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
	//string fn = generateData(36, 80);
	//readCostMatrix(fn);

	//tabu search for existing file
	readCostMatrix("dane36_50.txt");


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

	double best = solveTSP(10);

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
