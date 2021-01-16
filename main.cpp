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
#include <complex.h>
using namespace std;
int returnMinRow(vector<vector<double>> arr, int whichClass);
pair<int, int> returnMinMaxofColum(vector<vector<double>> arr, int whichcolumn);
int rankforme = -1;
void makeRelief(vector<vector<double>> &twodimensionalarr, int A, int P, int N, int iteration);
void makeManhantenDistance(vector<vector<double>> twodimensionalarr, int A, int iteration, vector<double> &manhattandistance, int nearestmissrow, int nearesthitrow, int M);
void printSlaveResult(vector<double> manhattandistance, int A, int T, set<int> &formasterp);
void make1Dto2D(vector<vector<double>> &twodimensionalarr2, vector<vector<double>> &twodimensionalarr, double local_dataset[], int A, int P, int N);
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

    int parts_size = (N / (P - 1)) * A + (N / (P - 1));
    double global_dataset[N * A + N + parts_size];
    double local_dataset[parts_size];
    for (int i = 0; i < parts_size; i++)
    {
        global_dataset[i] = 0;
    }
    if (rank == 0) // i am in master process
    {
        int index = parts_size;
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
                    global_dataset[index] = stod(var);
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
        global_dataset, parts_size, MPI_DOUBLE,
        local_dataset, parts_size, MPI_DOUBLE,
        0, MPI_COMM_WORLD);
    if (rank != 0)
    {

        set<int> formasterp;
        rankforme = rank;
        vector<int> selectedindexes;
        int nearestmissrow = 0, nearesthitrow = 0;
        vector<double> manhattandistance(A, 0);

        for (int iteration = 0; iteration < M; iteration++)
        {
            vector<vector<double>> twodimensionalarr;
            vector<vector<double>> twodimensionalarr2;
            make1Dto2D(twodimensionalarr2, twodimensionalarr, local_dataset, A, P, N);
            selectedindexes.push_back(iteration);
            makeRelief(twodimensionalarr2, A, P, N, iteration);
            nearesthitrow = returnMinRow(twodimensionalarr2, twodimensionalarr[iteration][twodimensionalarr[iteration].size() - 1]);
            nearestmissrow = returnMinRow(twodimensionalarr2, 1 - (twodimensionalarr[iteration][twodimensionalarr[iteration].size() - 1]));
            makeManhantenDistance(twodimensionalarr, A, iteration, manhattandistance, nearestmissrow, nearesthitrow, M);
        }
        printSlaveResult(manhattandistance, A, T, formasterp);

        int xx[1] = {5};
        MPI_Send(
            &xx,
            T,
            MPI_INT,
            0,
            0,
            MPI_COMM_WORLD);
    }
    if (rank == 0)
    {
        set<int> tt;
        int xx[1];
        int i=0;
        while(i < P)
        {
            MPI_Recv(
                &xx,
                T,
                MPI_INT,
                i,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
            cout << xx[0] << endl;
            for (auto m : tt)
                tt.insert(m);
            // }
            cout << "Master P0 : ";
            for (auto m : tt)
                cout << m << " ";
            cout << endl;
            i++;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }
}
void make1Dto2D(vector<vector<double>> &twodimensionalarr2, vector<vector<double>> &twodimensionalarr, double local_dataset[], int A, int P, int N)
{
    int yy = 0;
    for (int i = 0; i < (N / (P - 1)); i++) //converting one dimensional array into two dimensional array
    {
        vector<double> tempppp;
        for (int column = 0; column <= A; column++)
        {

            tempppp.push_back(local_dataset[yy]);
            yy++;
        }
        twodimensionalarr2.push_back(tempppp);
        twodimensionalarr.push_back(tempppp);
    }
}
void printSlaveResult(vector<double> manhattandistance, int A, int T, set<int> &formasterp)
{
    priority_queue<double> arr;
    vector<int> a;
    for (auto i : manhattandistance)
    {
        arr.push(i);
    }
    for (int i = 0; i < T; i++)
    {
        for (int j = 0; j < A; j++)
        {
            if (manhattandistance[j] == arr.top())
            {
                arr.pop();
                a.push_back(j);
                break;
            }
        }
    }
    cout << "Slave P" << rankforme << " :";
    sort(a.begin(), a.end());
    for (int i = 0; i < T; i++)
    {
        cout << a[i] << " ";
        formasterp.insert(a[i]);
    }
    cout << endl;
}
void makeManhantenDistance(vector<vector<double>> twodimensionalarr, int A, int iteration, vector<double> &manhattandistance, int nearestmissrow, int nearesthitrow, int M)
{
    for (int stepfour = 0; stepfour < A; stepfour++)
    {
        /*pair<minindex,maxindex>*/ pair<int, int> minmaxrow = returnMinMaxofColum(twodimensionalarr, stepfour);
        double temp = cabs(twodimensionalarr[minmaxrow.first][stepfour] - twodimensionalarr[minmaxrow.second][stepfour]) / M;
        double temp2 = cabs(twodimensionalarr[iteration][stepfour] - twodimensionalarr[nearestmissrow][stepfour]);
        double temp3 = cabs(twodimensionalarr[iteration][stepfour] - twodimensionalarr[nearesthitrow][stepfour]);
        manhattandistance[stepfour] += (temp2 / temp) - (temp3 / temp);
    }
}
void makeRelief(vector<vector<double>> &twodimensionalarr, int A, int P, int N, int iteration)
{
    for (int i = 0; i < (N / (P - 1)); i++)
    {
        for (int column = 0; column < A && i != iteration; column++)
        {
            twodimensionalarr[i][column] = cabs(twodimensionalarr[i][column] - twodimensionalarr[iteration][column]);
        }
    }
}
int returnMinRow(vector<vector<double>> arr, int whichClass)
{
    int minrow = INT16_MAX;
    int index = -1;
    int row = 0;
    for (row = 0; row < arr.size(); row++)
    {
        if (arr[row][arr[row].size() - 1] != whichClass)
            continue;
        int temp = 0;
        for (int column = 0; column < (arr[row].size() - 1); column++)
        {
            temp += arr[row][column];
        }
        if (temp < minrow)
        {
            index = row;
            minrow = temp;
        }
    }
    return index;
}

pair<int, int> returnMinMaxofColum(vector<vector<double>> arr, int whichcolumn)
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