#include <stdio.h>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <sstream>
using namespace std;
int main(int argc, char *argv[])
{

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    fstream input;
    int N, A, M, T;
    vector<vector<float>> dataset;
    pair<int, int> part_for_me [N];

    if (rank == 0) // i am in master process
    {
        string line, numberofprocess, n;
        input.open(argv[1]);
        getline(input, numberofprocess);
        stringstream oo(numberofprocess);
        oo >> n;
        //size = stoi(n);
        getline(input, line);
        stringstream ss(line);
        ss >> N >> A >> M >> T;
        int row;
        while (getline(input, line))
        {
            stringstream ss(line);
            do
            {
                if (ss.eof())
                    break;
                string temp;
                ss >> temp;
                try
                {
                    dataset[row].push_back(stof(temp));
                }
                catch (exception e)
                {
                    break;
                }
            } while (ss);
            row++;
        }
        int part = N / (size - 1);
        int index=0;
        for (int u = 0; u < (N - part); u += part)
        {
            part_for_me[index] = { u, u + part };
            index++;
        }
        for (int k = 0; k < size - 1; k++)
        {
        }
    }
    else
    {
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}