#include <iostream>
#include <stdlib.h>
//#include "tsptabusolver.h"
#include <ctime>
#include <string>
#include <fstream>
#include <limits>
using namespace std;


#define TABU_LENGTH 30
#define NUM_INTERATION 3000
#define PENAL_LONG_TERM 10
#define LONG_TERM_LENGTH 100
#define TIME_TRY 500

int N; // number of cities (vertices)
double **costMatrix; // cost of each edge
int **tabu_list, **tabu_f_list;
double score, bestSolverScore;
int *v, *foundSolution; // vertices, found path


void readCostMatrix(string fileName){
    ifstream file(fileName.c_str());

    if(!file.is_open()){
        perror ("The following error occurred");
        exit(0);
    }

    if(file.is_open()){
    	file >> N;

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

double getScore(int *v){ // cost of given path
	int score = 0;
	for(int i = 0; i < (N - 1); ++i){
		score += costMatrix[v[i]][v[i+1]];
	}
	return score += costMatrix[v[N-1]][v[0]];
}

void resetTabuList(){ // empty tabu list
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < N; ++j){
			tabu_list[i][j] = 0;
			tabu_f_list[i][j] = 0;
		}
	}
}

void initTabuList(){ // dynamically allocate tabu list
	tabu_list = new int*[N];
	tabu_f_list = new int*[N];
	for(int i = 0; i < N; ++i){
		tabu_f_list[i] = new int[N];
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


int* getBestNearbySolution(int* v, int it){ // search for neighbours
	int *v_temp = new int[N]; // copy of current solution

	for(int i = 0; i < N; ++i){
		v_temp[i]=v[i];
	}

	double bestScore = std::numeric_limits<double>::max();
	int vertexA = 0;
	int vertexB = 1;
	for(int i = 0; i < N; ++i){
		for(int j = (i+2); j < N; ++j){ // zmiana z j=(i+1)
			//swap for new solution
			swap(v_temp[i], v_temp[j]); //s->swapSolve(i,j);

			double currentScore = getScore(v_temp); //s->getScore();
			double penalScore = currentScore + PENAL_LONG_TERM * tabu_f_list[i][j];
			if( (bestScore > penalScore && tabu_list[i][j] <= it) || currentScore < bestSolverScore){
				vertexA = i;
				vertexB = j;
				bestScore = penalScore;
				tabu_list[i][j] = (it + TABU_LENGTH);
				tabu_list[j][i] = (it + TABU_LENGTH);
			}
			// back to orginal solution
			swap(v_temp[j], v_temp[i]);//s->swapSolve(j,i);
			if(tabu_f_list[i][j] > 0 && it > LONG_TERM_LENGTH) tabu_f_list[i][j] -= 1;
		}
	}
	tabu_f_list[vertexA][vertexB] += 2;
	swap(v_temp[vertexA], v_temp[vertexB]);//s->swapSolve( vertexA, vertexB );
	return v_temp;//s;
}

double solveTSP(int numCandidate){
	//Solution bestSolution(map);
	int *v_temp = new int[N];
	double bestSolutionScore = getScore(v);

	for(int loopCount = 0; loopCount < numCandidate; ++loopCount){
		initSolution();
		resetTabuList();
		//cout << "Init Score : " << s->getScore() << endl;
		int countTime = 0;
		bestSolverScore = std::numeric_limits<double>::max();

		for(int i = 0; i < NUM_INTERATION; ++i){
			v_temp = getBestNearbySolution(v, i);
			double score = getScore(v_temp);
			if(score < bestSolverScore){
				bestSolverScore = score;
				countTime = 0;

				if(bestSolverScore < bestSolutionScore){
					for(int j = 0; j < N; ++j){
						v[j]=v_temp[j];//bestSolution.set(j,s->getV(j));
						foundSolution[j] = v[j];
						//cout<<v[j]<<" ";
					}
					//cout<<endl;
					bestSolutionScore = bestSolverScore;
				}
			}else{
				++countTime;
				if(countTime > TIME_TRY){
					break;
				}
			}
		}

	}
	//cout << "Best score : " << bestSolutionScore << endl;
	//bestSolution.printPath();
	return bestSolutionScore;
}



int main(int argc, char* argv[]){
/*
	time_t start, end;

	time(&start);
	const clock_t begin_time = clock();

	TSPTabuSolver solver2("tsp0.txt");
	solver2.solve(6);	

	time(&end);
	cout << "Czas wykonania: " << difftime(end,start)<<endl;
	cout << "Czas procesora: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC <<endl;
*/
	/*
	TSPTabuSolver solver1("tsp1.txt");
	solver1.solve(5);				

				

	TSPTabuSolver solver3("tsp2.txt");
	solver3.solve(7);				
	*/

	readCostMatrix("file.txt");
	initSolution();
	initTabuList();

	int best = solveTSP(1);

	cout<<"Koñcowe rozwi¹zanie: ";
	for(int i = 0; i < N; ++i){
		cout<<foundSolution[i]<<" ";
	}
	cout<<endl<<"Koszt trasy: "<<best<<endl;
	return 0;
}
