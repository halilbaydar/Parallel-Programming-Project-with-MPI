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
void printSlaveResult(vector<double> manhattandistance, int A, int T, vector<int> &formasterp);
void make1Dto2D(vector<vector<double>> &twodimensionalarr2, vector<vector<double>> &twodimensionalarr, double local_dataset[], int A, int P, int N);
/*
In the main function all processors read first two line in order to determine the size of local data and global data arrays. 
Then Master process reads all N lines data and scatters it to all slaves by dividing the data equal size.
In the main function, if rank is not equal to rank 0 ,which means I am in slave processors, iterations steps are applied by 
slave processors.And finally results are written to console and also sent to master process. In order to wait all processors’ 
process I used MPI_Barrier(MPI_COMM_WORLD) to write all outputs to console before the master processor does. 
In master processors master processor receives all results from slave processors and writes them to console in ascending order. 
And all of them is terminated.
*/
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
    oo >> thenumberofprocess; //reads the number of processors
    int P = stoi(thenumberofprocess);
    int N, A, M, T;
    getline(input, line);
    stringstream ss(line);
    ss >> N >> A >> M >> T; //reads second line

    int parts_size = (N / (P - 1)) * A + (N / (P - 1));
    double global_dataset[N * A + N + parts_size];
    double local_dataset[parts_size];
    for (int i = 0; i < parts_size; i++)
    {
        global_dataset[i] = 0;
    }
    int ttemp[T];
    if (rank == 0) // i am in master process
    {
        //reads N data lines
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
    //scatter the all data to slave processors
    MPI_Scatter(
        global_dataset, parts_size, MPI_DOUBLE,
        local_dataset, parts_size, MPI_DOUBLE,
        0, MPI_COMM_WORLD);
    if (rank != 0)
    {
        //if rank is not equal to master processor, make iteration in given data
        vector<int> formasterp;
        rankforme = rank;
        vector<int> selectedindexes;
        int nearestmissrow = 0, nearesthitrow = 0;
        vector<double> manhattandistance(A, 0);

        //iteration for-loop
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
        //print result
        printSlaveResult(manhattandistance, A, T, formasterp);

        for (int uu = 0; uu < T; uu++)
            ttemp[uu] = formasterp[uu];
        //send T result to master processor    
        MPI_Send(
            &ttemp,
            T,
            MPI_INT,
            0,
            0,
            MPI_COMM_WORLD);
    }
    //wait all slave processors to continue properly
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        set<int> tt;
        int i = 1;
        cout << "Master P0 : ";
        while (i < P)
        {
            //receive all results from slave processors and print them in ascending order
            MPI_Recv(
                &ttemp,
                T,
                MPI_INT,
                i,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
            for (auto m : ttemp)
                tt.insert(m);
            // }

            i++;
        }
        for (auto m : tt)
            cout << m << " ";
        cout << endl;
    }
    //finilize
    MPI_Finalize();
    return 0;
}
//I broadcasts N lines data in one dimensional array and due to that I need to make 1D array to 2D array in order not to confuse me in processing all.
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
//This function writes the result to console and it is called only in slave processes.
void printSlaveResult(vector<double> manhattandistance, int A, int T, vector<int> &formasterp)
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
        formasterp.push_back(a[i]);
    }
    cout << endl;
}
//This function applies the manhanten distance formula with given parameters
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
/*
This function makes absolute subtract operation between the selected instance and  
all other instance expext itself to help in finding the nearest hit and the nearest miss.
*/
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
/*
This function iterates the given data according to considered class type and sums all 
features in one class and compares other sums of same class. Finally returns given class 
row of the minimum sum. That is, If class parameter is same with the iterated class, the 
return value is the nearest hit, otherwise the nearest miss.*/
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

/*
This function returns the indicies of the max and min value in given column to apply manhanten distance algowithm. */
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