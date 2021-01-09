#include <stdio.h>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <vector>
#include <queue>
#include <map>
#include <bits/stdc++.h>
using namespace std;
int returnMinRow(vector<vector<float>> arr, int whichClass);
pair<int, int> returnMinMaxofColum(vector<vector<float>> arr, int whichcolumn);
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

    if (rank != 0)
    {

        MPI_Gather(
            local_dataset, parts_size, MPI_FLOAT,
            global_dataset, parts_size, MPI_FLOAT,
            1, MPI_COMM_WORLD);

        vector<int> selectedindexes;
        int nearestmissrow = 0, nearesthitrow = 0;
        float manhattandistance[A];
        for (int iteration = 0; iteration < M; iteration++)
        {
            vector<vector<float>> twodimensionalarr;
            vector<vector<float>> twodimensionalarr2;
            int yy = 0;
            for (int i = 0; i < (N / (P - 1)); i++) //converting one dimensional array into two dimensional array
            {
                for (int column = 0; column <= A; column++)
                {
                    twodimensionalarr[i][column] = local_dataset[yy];
                    twodimensionalarr2[i][column] = local_dataset[yy];
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
            for (int i = 0; i < (N / (P - 1)); i++)
            {
                for (int column = 0; column < A && i != randomrow; column++)
                {
                    twodimensionalarr[i][column] -= twodimensionalarr[randomrow][column];
                }
            }

            nearesthitrow = returnMinRow(twodimensionalarr, twodimensionalarr[randomrow][twodimensionalarr[randomrow].size() - 1]);
            nearestmissrow = returnMinRow(twodimensionalarr, 1 - (twodimensionalarr[randomrow][twodimensionalarr[randomrow].size() - 1]));

            for (int stepfour = 0; stepfour < A; stepfour++)
            {
                /*pair<minindex,maxindex>*/ pair<int, int> minmaxrow = returnMinMaxofColum(twodimensionalarr2, stepfour);
                float temp = abs(twodimensionalarr2[minmaxrow.first][stepfour] - twodimensionalarr[minmaxrow.second][stepfour]) / M;
                float temp2 = abs(twodimensionalarr2[randomrow][stepfour] - twodimensionalarr[nearestmissrow][stepfour]);
                float temp3 = abs(twodimensionalarr2[randomrow][stepfour] - twodimensionalarr[nearesthitrow][stepfour]);
                manhattandistance[stepfour] -= (temp3 / temp) + (temp2 / temp);
            }
        }
        queue<float> manhattandistance2;
        for (int p = 0; p < A; p++)
        {
            manhattandistance2.push(manhattandistance[p]);
        }
        MPI_Send(
            &manhattandistance2,
            1,
            MPI_FLOAT,
            0,
            0,
            MPI_COMM_WORLD);
        for (int k = 0; k < T; k++)
        {
            cout << manhattandistance2.front();
            manhattandistance2.pop();
        }
    }
    if (rank == 0)
    {
        set<float> results;
        for (int i = 1; i <= P; i++)
        {
            queue<float> receive;
            MPI_Recv(
                &receive,
                1,
                MPI_FLOAT,
                1,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
            for (int t = 0; t < T; t++)
            {
                results.insert(receive.front());
                receive.pop();
            }
        }
        for (set<float>::iterator t = results.begin(); t != results.end(); t++)
            cout << *t < endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
int returnMinRow(vector<vector<float>> arr, int whichClass)
{
    int minrow = INT16_MAX;
    int index = -1;
    for (int row = 0; row < arr.size(); row++)
    {
        int temp = 0;
        for (int column = 0; column < arr[row].size() - 1 && arr[row][arr[row].size() - 1] == whichClass; column++)
        {
            temp += arr[row][column];
        }
        if (temp < minrow && arr[row][arr[row].size() - 1] == whichClass)
        {
            index = row;
            minrow = temp;
        }
    }
    return index;
}

pair<int, int> returnMinMaxofColum(vector<vector<float>> arr, int whichcolumn)
{
    double Min = INT_MAX, Max = INT_MIN;
    int minindex, maxindex;
    for (int row = 0; row < arr.size(); row++)
    {
        if (Min > arr[row].at(whichcolumn))
        {
            Min = arr[row][whichcolumn];
            minindex = row;
        }
        if (Max < arr[row][whichcolumn])
        {
            Max = arr[row][whichcolumn];
            maxindex = row;
        }
    }
    return {minindex, maxindex};
}
