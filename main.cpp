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
int returnMinRow(vector<vector<double>> arr, int whichClass);
pair<int, int> returnMinMaxofColum(vector<vector<double>> arr, int whichcolumn);
int rankforme = -1;
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
                    //cout << "arr size küçük";
                    break;
                }
                index++;
            } while (ss);
        }

        ////cout<<"okuma başarılı"<<endl;
    }

    MPI_Scatter(
        global_dataset, parts_size, MPI_DOUBLE,
        local_dataset, parts_size, MPI_DOUBLE,
        0, MPI_COMM_WORLD);
    //cout<<"distribute başarılı"<<endl;
    if (rank != 0)
    {
        rankforme = rank;
        //cout<<rank<<" işlem yapıyor"<<endl;
        vector<int> selectedindexes;
        int nearestmissrow = 0, nearesthitrow = 0;
        double manhattandistance[A];
        for (int iteration = 0; iteration < M; iteration++)
        {
            //cout<<"iteration yapılıyor rank: "<<rank<<" iteration : "<<iteration<<endl;
            vector<vector<double>> twodimensionalarr;
            vector<vector<double>> twodimensionalarr2;
            int yy = 0;
            for (int i = 0; i < (N / (P - 1)); i++) //converting one dimensional array into two dimensional array
            {
                vector<double> tempppp;
                for (int column = 0; column <= A; column++)
                {
                    tempppp.push_back(local_dataset[yy]);
                    tempppp.push_back(local_dataset[yy]);
                    yy++;
                }
                twodimensionalarr2.push_back(tempppp);
                twodimensionalarr.push_back(tempppp);
            }
            //cout<<"random row seçiliyor"<<endl;
            int randomrow = rand() % (N / (P - 1)) + 1; //1 den verilen row sayısına kadar
            vector<int>::iterator myit = find(selectedindexes.begin(), selectedindexes.end(), randomrow);
            while (myit != selectedindexes.end())
            {
                randomrow = rand() % (N / (P - 1)) + 1; //1 den verilen row sayısına kadar
                myit = find(selectedindexes.begin(), selectedindexes.end(), randomrow);
            }
            selectedindexes.push_back(randomrow);
            //cout<<"rank :"<<rank <<"rowlar random rowdan çıkarılıyor"<<endl;
            for (int i = 0; i < (N / (P - 1)); i++)
            {
                for (int column = 0; column < A && i != randomrow; column++)
                {
                    twodimensionalarr[i][column] -= twodimensionalarr[randomrow][column];
                }
            }
            //cout<<"rank :"<<rank<<" narest hit ve miss hesaplanıyor"<<endl;
            nearesthitrow = returnMinRow(twodimensionalarr, twodimensionalarr[randomrow][twodimensionalarr[randomrow].size() - 1]);
            nearestmissrow = returnMinRow(twodimensionalarr, 1 - (twodimensionalarr[randomrow][twodimensionalarr[randomrow].size() - 1]));

            //cout<<"nearestmissrow :"<<nearestmissrow << " nearesthitrow :" <<nearesthitrow<<endl;
            for (int stepfour = 0; stepfour < A; stepfour++)
            {
                /*pair<minindex,maxindex>*/ pair<int, int> minmaxrow = returnMinMaxofColum(twodimensionalarr2, stepfour);
                //cout<<"minindex : "<<minmaxrow.first<<" maxindex :" <<minmaxrow.second<<endl;
                double temp = abs(twodimensionalarr2[minmaxrow.first][stepfour] - twodimensionalarr[minmaxrow.second][stepfour]) / M;
                double temp2 = abs(twodimensionalarr2[randomrow][stepfour] - twodimensionalarr[nearestmissrow][stepfour]);
                double temp3 = abs(twodimensionalarr2[randomrow][stepfour] - twodimensionalarr[nearesthitrow][stepfour]);
                manhattandistance[stepfour] = manhattandistance[stepfour]-(temp3 / temp) + (temp2 / temp);
            }
            //cout<<"rank :"<<rank<<" manhantan distance hesaplandı "<<endl;
        }
        for(auto i:manhattandistance)
            cout<<i<<" ";
        cout<<endl;
        set<double> manhattandistance2;
        for (int p = 0; p < A; p++)
        {
            manhattandistance2.insert(manhattandistance[p]);
        }
        cout << "Slave P" << rank << " :";
        int oo = 0;
        for (set<double>::iterator it = manhattandistance2.begin(); it != manhattandistance2.end(); it++)
        {
            for (int j = 0; j < A; j++)
            {
                if (*it == manhattandistance[j])
                {
                    cout << j << " ";
                    break;
                }
            }
            oo++;
            if (oo == T)
                break;
        }
        cout << endl;
    }
    /*
   if(rank==0){
       cout<<"rank"<<rank<<endl;
        for(int i=0; i<parts_size; i++)
            cout<<local_dataset[i]<<" ";
        cout<<endl;
        
       cout<<parts_size;
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
int returnMinRow(vector<vector<double>> arr, int whichClass)
{
    int minrow = INT16_MAX;
    int index = -1;
    int row = 0;
    for (row = 0; row < arr.size(); row++)
    {
        int temp = 0;
        for (int column = 0; column < (arr[row].size() - 1) && arr[row][arr[row].size() - 1] == whichClass; column++)
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
