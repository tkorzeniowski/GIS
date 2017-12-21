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
double score=0.0, bestSolverScore=0.0;
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


//int* getBestNearbySolution(int* v, int it){ // search for neighbours
int* getBestNearbySolution(int it){ // search for neighbours
	int *v_temp = new int[N]; // copy of current solution

	for(int i = 0; i < N; ++i){
		v_temp[i]=v[i];
	}

	double bestScore = std::numeric_limits<double>::max();
	int vertexA = 0;
	int vertexB = 1;
	for(int i = 0; i < N; ++i){
		for(int j = (i+1); j < N; ++j){ // zmiana z j=(i+1)
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
			//v_temp = getBestNearbySolution(v, i);
			v_temp = getBestNearbySolution(i);


//			for(int h = 0; h < N; ++h){
//					//cout<<v_temp[h]<<" ";
//				}
//			//cout<<endl;


			double score = getScore(v_temp);
//			cout<<i<<" score "<<score<< " bestSolverScore "<<bestSolverScore<<" bestSolScore "<<bestSolutionScore<<endl;

//			double w=0;
//			for(int g = 0; g < (N - 1); ++g){
//						w += costMatrix[v_temp[g]][v_temp[g+1]];
//						//cout<<w<< " " <<costMatrix[v[g]][v[g+1]]<<endl;
//					}
//				w += costMatrix[v_temp[N-1]][v_temp[0]];
//			cout<<"koszt obliczony "<<w<<endl;


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

	readCostMatrix("test36.txt");
	initSolution();
//	double wynik = 0;
//	for(int i = 0; i < (N - 1); ++i){
//				wynik += costMatrix[v[i]][v[i+1]];
//				cout<<wynik<< " " <<costMatrix[v[i]][v[i+1]]<<endl;
//			}
//		wynik += costMatrix[v[N-1]][v[0]];
	initTabuList();

	double best = solveTSP(2);

	cout<<"Koñcowe rozwi¹zanie: ";
	for(int i = 0; i < N; ++i){
		cout<<foundSolution[i]<<" ";
	}
	cout<<endl<<"Koszt trasy: "<<best<<endl;

//	wynik=0;
//	for(int i = 0; i < (N - 1); ++i){
//			wynik += costMatrix[foundSolution[i]][foundSolution[i+1]];
//			cout<<wynik<< " " <<costMatrix[foundSolution[i]][foundSolution[i+1]]<<endl;
//		}
//	wynik += costMatrix[foundSolution[N-1]][foundSolution[0]];
//	cout<<wynik<< " " <<costMatrix[foundSolution[N-1]][foundSolution[0]]<<endl;
//	cout<<"wynik = "<<wynik<<endl;

	return 0;
}
