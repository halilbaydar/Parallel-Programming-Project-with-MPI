#include <stdio.h>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <vector>
#include <queue>
#include <map>
using namespace std;
int main(int argc, char *argv[])
{

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    fstream input;
    string line, numberofprocess, thenumberofprocess;
    input.open(argv[1]);
    getline(input, numberofprocess);
    stringstream oo(numberofprocess);
    oo >> thenumberofprocess;
    int P = stoi(thenumberofprocess);
    int N, A, M, T;
    getline(input, line);
    stringstream ss(line);
    ss >> N >> A >> M >> T;

    float global_dataset[N * A + N];
    int parts_size = pow((N / (P - 1)), 2) + (N / (P - 1));
    float local_dataset[parts_size];
    if (rank == 0) // i am in master process
    {
        int index;
        while (getline(input, line))
        {
            stringstream ss(line);
            do
            {
                if (ss.eof())
                    break;
                string var;
                ss >> var;
                try
                {
                    global_dataset[index] = stoi(var);
                }
                catch (exception e)
                {
                    break;
                }
                index++;
            } while (ss);
        }
    }
    MPI_Scatter(
        global_dataset, parts_size, MPI_FLOAT,
        local_dataset, parts_size, MPI_FLOAT,
        1, MPI_COMM_WORLD);

    MPI_Gather(
        local_dataset, parts_size, MPI_FLOAT,
        global_dataset, parts_size, MPI_FLOAT,
        1, MPI_COMM_WORLD);

    if (rank != 0)
    {
        vector<int> selectedindexes;
        map<float,int> resultmap;
        int nearestmissrow=0,nearesthitrow=0;
        for (int iteration = 0; iteration < M; iteration++)
        {
            float twodimensionalarr[N / (P - 1))][A];
            int yy = 0;
            for (int i = 0; i < (N / (P - 1)); i++) //converting one dimensional array into two dimensional array
            {
                for (int column = 0; column <= A; column++)
                {
                    twodimensionalarr[i][column] = local_dataset[yy];
                    yy++;
                }
            }
            int randomrow = rand() % (N / (P - 1)) + 1; //1 den verilen row sayısına kadar
            vector<int>::iterator myit = find(selectedindexes.begin(), selectedindexes.end(), randomrow);
            while (myit != selectedindexes.end())
            {
                randomrow = rand() % (N / (P - 1)) + 1; //1 den verilen row sayısına kadar
                myit = find(selectedindexes.begin(), selectedindexes.end(), randomrow);
            }
            selectedindexes.push_back(randomrow);
            for (int i = 0; i < (N / (P - 1)); i++) //converting one dimensional array into two dimensional array
            {
                for (int column = 0; column < A && i != randomrow; column++)
                {
                    twodimensionalarr[i][column] -= twodimensionalarr[randomrow][column];
                }
            }
            for (int i = 0; i < (N / (P - 1)); i++) //converting one dimensional array into two dimensional array
            {
                int thesumofrow=0;
                for (int column = 0; column < A && i != randomrow; column++)
                {
                    thesumofrow= twodimensionalarr[i][column];
                }
                resultmap[thesumofrow]=i;
            }
            //for nearesetmissrow
            int indextemp=0;
            for(map<float,int>::iterator item=resultmap.begin(); item!=resultmap.end(); item++){
                    if(item->second==twodimensionalarr[indextemp][A]){
                        nearestmissrow=indextemp;
                        break;
                    }
            }
            //for nearesthitrow
            int indextemp=0;
            for(map<float,int>::iterator item=resultmap.begin(); item!=resultmap.end(); item++){
                    if(item->second!=twodimensionalarr[indextemp][A]){
                        nearesthitrow=indextemp;
                        break;
                    }
            }
    
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

/*for (int i = 0; i < (N / (P - 1)); i++)
            {
                for (int column = 0; column <= A && (A * i + i != ((randomrow - 1) * A + (randomrow - 1))); column++)
                {
                    new_local_dataset[A * i + i + (i == 0 ? column : (column + 1))] =
                        abs(new_local_dataset[A * i + i + (i == 0 ? column : (column + 1))] -
                            local_dataset[((randomrow - 1) * A + (randomrow - 1)) + (column + 1)]) //randomly selected row
                }
            }*/