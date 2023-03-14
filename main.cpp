#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <ctime>
#include <chrono>
#include <limits>

using namespace std;

class Problem
{
public:
    // function to get the number of cities
    int getNumCities(string fileName)
    {
        int numCities = 0;
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

    // function to get the type of the problem
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

    // function to read the file and store the coordinates of the cities
    double **readFromFile(string fileName, double **cities)
    {
        int numCities = 0;
        // double **cities;
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
            cities = new double *[numCities];

            for (int i = 0; i < numCities; i++)
            {
                cities[i] = new double[2];
            }

            // if str is "NODE_COORD_SECTION"
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
    // function to generate a random solution
    int *getRandomSolution(int *arr, int size)
    {
        // Seed the random number generator
        srand(time(NULL));

        // Iterate through the array starting from the last element
        for (int i = size - 1; i > 0; i--)
        {
            // Generate a random index between 0 and i (inclusive)
            int j = rand() % (i + 1);

            // Swap the current element with a random element before it
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
        return arr;
    }

    // function to calculate euclidean distance between two points
    double **euclideanDistance(double **cities, double **distances, int numCities)
    {
        distances = new double *[numCities];
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
    double **geoDistance(double **cities, double **distances, int numCities)
    {
        distances = new double *[numCities];
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
    // function to calculate distance between two points
    double **calculateDistance(double **cities, double **distances, int numCities, string type)
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
    double getCost(int *sol, double **distances, int numCities)
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
    // random search algorithm
    int *randomAlgorithm(int *best_solution, Solution solution, double **distances, int numCities, double min_cost, time_t start_time, int time_limit)
    {
        int iter_number = 0;
        do
        {
            iter_number++;
            int *new_solution = solution.getRandomSolution(best_solution, numCities);
            double new_cost = solution.getCost(new_solution, distances, numCities);
            if (new_cost < min_cost)
            {
                min_cost = new_cost;
                best_solution = new_solution;
            }
        } while ((time(NULL) - start_time) < time_limit);

        cout << "Number of iterations: " << iter_number << endl;

        return best_solution;
    }
    // random walk algorithm
    int *randomWalkAlgorithm(int *best_solution, Solution solution, double **distances, int numCities, double min_cost, time_t start_time, int time_limit)
    {
        int iter_number = 0;
        do
        {
            iter_number++;
            // shuffle random part of the solution
            int start = rand() % numCities;
            int end = start + rand() % (numCities - start);

            // copy the best solution to a new solution
            int *new_solution = new int[numCities];
            for (int i = 0; i < numCities; i++)
            {
                new_solution[i] = best_solution[i];
            }

            // shuffle the array between start and end
            for (int i = start; i < end; i++)
            {
                int j = rand() % (end - start) + start;
                int temp = new_solution[i];
                new_solution[i] = new_solution[j];
                new_solution[j] = temp;
            }

            double new_cost = solution.getCost(new_solution, distances, numCities);
            if (new_cost < min_cost)
            {
                min_cost = new_cost;
                best_solution = new_solution;
            }
        } while ((time(NULL) - start_time) < time_limit);

        cout << "Number of iterations:" << iter_number << endl;

        return best_solution;
    }

    // nearest neighbor algorithm
    int *nearestNeighborAlgorithm(int *best_solution, double **distances, int numCities)
    {
        auto start_time = chrono::high_resolution_clock::now();
        int *new_solution = new int[numCities];
        new_solution[0] = best_solution[0];
        int *visited = new int[numCities];
        for (int i = 0; i < numCities; i++)
        {
            visited[i] = 0;
        }
        visited[best_solution[0]] = 1;
        for (int i = 1; i < numCities; i++)
        {
            double min_dist = numeric_limits<double>::infinity();;
            int min_index = 0;
            for (int j = 0; j < numCities; j++)
            {
                if (visited[j] == 0 && distances[new_solution[i - 1]][j] < min_dist)
                {
                    min_dist = distances[new_solution[i - 1]][j];
                    min_index = j;
                }
            }
            new_solution[i] = min_index;
            visited[min_index] = 1;
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "Time taken by function: " << chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << endl;
        return new_solution;
    }
};

int main()
{
    // read data from file
    string fileName = "tsp_instances/kroA100.tsp"; // state file name
    Problem problem;
    double **cities;
    double **distances;
    int numCities = problem.getNumCities(fileName); // get number of cities
    string type = problem.getType(fileName);        // get type of distance
    cout << "Number of cities:" << numCities << endl;
    cout << "Type of distance:" << type << endl;
    cities = problem.readFromFile(fileName, cities); // read cities from file

    // -----------------------------------------

    // calculate all distances, get initial solution and its cost
    Solution solution;
    distances = solution.calculateDistance(cities, distances, numCities, type); // calculate distance between cities

    int *sol = new int[numCities];
    for (int i = 0; i < numCities; i++)
    {
        sol[i] = i;
    }
    sol = solution.getRandomSolution(sol, numCities);          // get initial solution
    double cost = solution.getCost(sol, distances, numCities); // get cost of initial solution
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------

    // perform different algorithms
    Algorithm algorithm;
    // RANDOM
    cout << endl
         << "RANDOM" << endl;
    double min_cost = cost; // cost of initial solution
    int *best_solution;
    int time_limit = 3; // time limit in seconds
    time_t start_time = time(NULL);
    best_solution = algorithm.randomAlgorithm(sol, solution, distances, numCities, min_cost, start_time, time_limit);
    // print best solution
    for (int i = 0; i < numCities; i++)
    {
        cout << best_solution[i] << " ";
    }
    cout << endl;
    cout << solution.getCost(best_solution, distances, numCities) << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // RANDOM WALK
    cout << endl
         << "RANDOM WALK" << endl;
    min_cost = cost; // cost of initial solution
    time_limit = 3;  // time limit in seconds
    start_time = time(NULL);
    best_solution = algorithm.randomWalkAlgorithm(sol, solution, distances, numCities, min_cost, start_time, time_limit);
    // print best solution
    for (int i = 0; i < numCities; i++)
    {
        cout << best_solution[i] << " ";
    }
    cout << endl;
    cout << solution.getCost(best_solution, distances, numCities) << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // NEAREST NEIGHBOR
    cout << endl
         << "NEAREST NEIGHBOR" << endl;
    best_solution = algorithm.nearestNeighborAlgorithm(sol, distances, numCities);
    // print best solution
    for (int i = 0; i < numCities; i++)
    {
        cout << best_solution[i] << " ";
    }
    cout << endl;
    cout << solution.getCost(best_solution, distances, numCities) << endl;
    cout << "-----------------------------------------" << endl;

    // free the dynamically allocated memory
    for (int i = 0; i <= numCities; i++)
    {
        delete[] cities[i];
        delete[] distances[i];
    }
    delete[] cities;
    delete[] distances;
    delete[] sol;
    delete[] best_solution;

    return 0;
}
