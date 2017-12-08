#include <iostream>
#include <string>
#include <fstream>
using namespace std;

void readCostM()
{
    ifstream file("file.txt");
    if(file.is_open())
    {
    	int N = 0;
    	file >> N;
        double myArray[N][N];

        for(int i = 0; i < N; ++i)
        {
        	for(int j=0; j<N; ++j)
        	{
        		file >> myArray[i][j];
        	}

        }

        for (int i=0; i<N; ++i){
        	cout<<myArray[i][i]<<endl;
        }
    }



}
