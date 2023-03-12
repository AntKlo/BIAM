#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <ctime>

using namespace std;


class Problem
{
public:
    int getNumCities(string fileName)
    {
        int numCities = 0;
        //double **cities;
        ifstream file(fileName);
        string str;

        while (getline(file, str))
        {
            // if str starts with "DIMENSIONS"
            if (str.rfind("DIMENSION", 0) == 0)
            {
                // get the number of cities
                int pos = str.find(":");
                string num = str.substr(pos + 1);
                numCities = stoi(num);
                return numCities;
            }
        }
    }

    string getType(string fileName)
    {
        string type;
        ifstream file(fileName);
        string str;

        while (getline(file, str))
        {
            // if str starts with "DIMENSIONS"
            if (str.rfind("EDGE_WEIGHT_TYPE", 0) == 0)
            {
                // get the number of cities
                int pos = str.find(":");
                type = str.substr(pos + 2);
                return type;
            }
        }
    }

    double** readFromFile(string fileName, double** cities)
    {
        int numCities = 0;
        //double **cities;
        ifstream file(fileName);
        string str;

        while (getline(file, str))
        {
            // if str starts with "DIMENSIONS"
            if (str.rfind("DIMENSION", 0) == 0)
            {
                // get the number of cities
                int pos = str.find(":");
                string num = str.substr(pos + 1);
                numCities = stoi(num);
            }

            // create 2d array to store the coordinates of the cities
            cities = new double*[numCities]; // dynamically allocate memory for the array

            for(int i=0; i<numCities; i++)
            {
                cities[i] = new double[2];
            }

            // str is "NODE_COORD_SECTION"
            if (str.rfind("NODE_COORD_SECTION", 0) == 0)
            {
                // read the coordinates of the cities
                for (int i = 0; i < numCities; i++)
                {
                    cities[i] = new double[2];
                    getline(file, str);
                    if (isspace(str[0]))
                    {
                        str = str.substr(1);
                    }
                    int pos = str.find(" ");
                    string str_copy = str.substr(pos + 1);
                    int pos_end = str_copy.find(" ");
                    string x = str_copy.substr(0, pos_end);
                    pos = str_copy.find(" ");
                    str_copy = str_copy.substr(pos + 1);
                    string y = str_copy.substr(0);
                    cities[i][0] = stod(x);
                    cities[i][1] = stod(y);
                }
                return cities;
            }
        }
    }
};


class Solution
{
public:
    // function to randomize vector
    int* getInitialSolution(int* arr, int n)
    {
        arr = new int[n];
        for (int i = 0; i < n; i++)
        {
            arr[i] = i;
        }
        // Start from the last element and swap
        // one by one. We don't need to run for
        // the first element that's why i > 0
        for (int i = n - 1; i > 0; i--)
        {
            // Pick a random index from 0 to i
            int j = rand() % (i + 1);

            // Swap arr[i] with the element
            // at random index
            arr[j] = arr[i] + arr[j];
            arr[i] = arr[j] - arr[i];
            arr[j] = arr[j] - arr[i];
        }
        return arr;
    }

    // function to calculate euclidean distance between two points
    double** euclideanDistance(double** cities, double** distances, int numCities)
    {
        distances = new double*[numCities];
        for (int i = 0; i < numCities; i++)
        {
            distances[i] = new double[numCities];
            for (int j = 0; j < numCities; j++)
            {
                double x1 = cities[i][0];
                double y1 = cities[i][1];
                double x2 = cities[j][0];
                double y2 = cities[j][1];
                double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
                distances[i][j] = distance;
            }
        }
        return distances;
    }

    // function to calculate geo distance between two points
    double** geoDistance(double** cities, double** distances, int numCities)
    {
        distances = new double*[numCities];
        for (int i = 0; i < numCities; i++)
        {
            distances[i] = new double[numCities];
            for (int j = 0; j < numCities; j++)
            {
                double x1 = cities[i][0];
                double y1 = cities[i][1];
                double x2 = cities[j][0];
                double y2 = cities[j][1];
                double distance = acos(sin(x1) * sin(x2) + cos(x1) * cos(x2) * cos(y2 - y1)) * 6378.388;
                distances[i][j] = distance;
            }
        }
        return distances;
    }
    // function that calculates euclidean distance between two points
    double** calculateDistance(double** cities, double** distances, int numCities, string type)
    {
        if (type == "EUC_2D")
        {
            distances = euclideanDistance(cities, distances, numCities);
        }
        else if (type == "GEO")
        {
            distances = geoDistance(cities, distances, numCities);
        }
        else
        {
            cout << "Invalid type" << endl;
        }
        return distances;
    }

    // function to calculate the cost of a solution
    double getCost(int* sol, double** distances, int numCities)
    {
        double cost = 0;
        for (int i = 0; i < numCities - 1; i++)
        {
            cost += distances[sol[i]][sol[i + 1]];
        }
        cost += distances[sol[numCities - 1]][sol[0]];
        return cost;
    }
};


class Algorithm
{
public:
    int* getRandomSolution(int* best_solution, int n, Solution solution, double** distances, int numCities, double min_cost, time_t start_time, int time_limit)
    {

        int* new_solution = solution.getInitialSolution(new_solution, n);
        double new_cost = solution.getCost(new_solution, distances, numCities);

        while ((time(0) - start_time) < time_limit)
        {
            // if the new solution is better than the current solution
            if (new_cost < min_cost)
            {
                int* best_solution = new_solution;
                cout << "New best solution found: " << new_cost << endl;
                return getRandomSolution(best_solution, n, solution, distances, numCities, new_cost, start_time, time_limit);
            }
            else
            {
                return getRandomSolution(best_solution, n, solution, distances, numCities, min_cost, start_time, time_limit);
            }
        }
        return best_solution;
    }    
};


int main()
{
    string fileName = "tsp_instances/kroA100.tsp";
    Problem problem;
    double** cities;
    double** distances;
    int numCities = problem.getNumCities(fileName);
    string type = problem.getType(fileName);
    cout << numCities << endl;
    cout << type << endl;
    cities = problem.readFromFile(fileName, cities);

    Solution solution;
    distances = solution.calculateDistance(cities, distances, numCities, type);
    
    
    int *sol;
    
    sol = solution.getInitialSolution(sol, numCities);
    double cost = solution.getCost(sol, distances, numCities);

    for (int i = 0; i < numCities; i++)
    {
        cout << sol[i] << " ";
    }

    cout << endl;
    cout << cost << endl;

    cout << "Random Solution" << endl;
    Algorithm algorithm;
    time_t start_time = time(0);
    int time_limit = 3;
    double min_cost = cost;
    int* best_solution;
    best_solution = algorithm.getRandomSolution(sol, numCities, solution, distances, numCities, min_cost, start_time, time_limit);
    // print best solution
    for (int i = 0; i < numCities; i++)
    {
        cout << best_solution[i] << " ";
    }
    cout << endl;
    cout << solution.getCost(best_solution, distances, numCities) << endl;


    // free the dynamically allocated memory
    for(int i=0; i<=numCities; i++)
    {
        delete[] cities[i];
        delete[] distances[i];
    }
    delete[] cities;
    delete[] distances;
    delete[] sol;

    return 0;
}