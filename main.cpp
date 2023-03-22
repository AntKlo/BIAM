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

void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

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

    // function to save 10 solutions from array and their costs to a file
    void saveToFile(string fileName, int **solutions, double *costs, int numCities)
    {
        ofstream file;
        file.open(fileName);
        for (int i = 0; i < 10; i++)
        {
            for (int j = 0; j < numCities; j++)
            {
                file << solutions[i][j] << " ";
            }
            file << endl;
            file << costs[i] << endl;
        }
        file.close();
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

        cout << "Number of iterations: " << iter_number << endl;

        return best_solution;
    }

    // nearest neighbor algorithm
    int *nearestNeighborAlgorithm(double **distances, int numCities)
    {
        auto start_time = chrono::high_resolution_clock::now();
        int *new_solution = new int[numCities];
        new_solution[0] = rand() % numCities;
        bool *visited = new bool[numCities];
        for (int i = 0; i < numCities; i++)
        {
            visited[i] = false;
        }
        visited[new_solution[0]] = true;
        for (int i = 1; i < numCities; i++)
        {
            double min_dist = numeric_limits<double>::infinity();
            int min_index = 0;
            for (int j = 0; j < numCities; j++)
            {
                if (visited[j] == false && distances[new_solution[i - 1]][j] < min_dist)
                {
                    min_dist = distances[new_solution[i - 1]][j];
                    min_index = j;
                }
            }
            new_solution[i] = min_index;
            visited[min_index] = true;
        }
        auto end_time = chrono::high_resolution_clock::now();
        // cout << "Time taken by function: " << chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << endl;
        delete[] visited;
        return new_solution;
    }

    // greedy algorithm
    int *greedyAlgorithm(int *initial_solution, double **distances, int numCities)
    {
        Solution solution;
        auto start_time = chrono::high_resolution_clock::now();
        // local search
        int *new_solution = new int[numCities];
        for (int i = 0; i < numCities; i++)
        {
            new_solution[i] = initial_solution[i];
        }

        // greedy improvement
        bool improved = true;
        int iter_number = 0;
        int solutions_visited = 0;
        while (improved)
        {
            iter_number++;
            improved = false;
            for (int i = 1; i < numCities - 1; i++)
            {
                for (int k = i + 1; k < numCities - 1; k++)
                {
                    solutions_visited++;
                    double delta = distances[new_solution[i - 1]][new_solution[k]] + distances[new_solution[i]][new_solution[k + 1]] - distances[new_solution[i - 1]][new_solution[i]] - distances[new_solution[k]][new_solution[k + 1]];
                    if (delta < 0)
                    {
                        int *temp_solution = new int[numCities];
                        for (int i = 0; i < numCities; i++)
                        {
                            temp_solution[i] = new_solution[i];
                        }

                        // two_opt_swap(i, k);
                        int m = i;
                        int n = k;
                        while (m < n)
                        {
                            int temp = temp_solution[m];
                            temp_solution[m] = temp_solution[n];
                            temp_solution[n] = temp;
                            m++;
                            n--;
                        }
                        if (solution.getCost(temp_solution, distances, numCities) < solution.getCost(new_solution, distances, numCities))
                        {
                            for (int i = 0; i < numCities; i++)
                            {
                                new_solution[i] = temp_solution[i];
                            }
                            improved = true;
                        }
                        if (improved)
                        {
                            break;
                        }
                    }
                }
                if (improved)
                {
                    break;
                }
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        // cout << "Time taken by function: " << chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << endl;
        // cout << "Number of iterations: " << iter_number << endl;
        // cout << "Number of solutions visited: " << solutions_visited << endl;
        return new_solution;
    }

    // steepest algorithm
    int *steepestAlgorithm(int *initial_solution, double **distances, int numCities)
    {
        Solution solution;
        auto start_time = chrono::high_resolution_clock::now();
        // local search
        int *best_solution = new int[numCities];
        for (int i = 0; i < numCities; i++)
        {
            best_solution[i] = initial_solution[i];
        }
        int *current_solution = new int[numCities];
        for (int i = 0; i < numCities; i++)
        {
            current_solution[i] = initial_solution[i];
        }

        bool improved = true;
        int iter_number = 0;
        int solutions_visited = 0;
        while (improved)
        {
            improved = false;
            iter_number++;
            double best_delta = 0;
            int best_i = -1;
            int best_k = -1;

            // iterate over all possible moves
            for (int i = 1; i < numCities - 1; i++)
            {
                for (int k = i + 1; k < numCities - 1; k++)
                {
                    solutions_visited++;
                    double delta = distances[current_solution[i - 1]][current_solution[k]] + distances[current_solution[i]][current_solution[k + 1]] - distances[current_solution[i - 1]][current_solution[i]] - distances[current_solution[k]][current_solution[k + 1]];
                    if (delta < best_delta)
                    {
                        best_delta = delta;
                        best_i = i;
                        best_k = k;
                    }
                }
            }

            // apply the best move if it improves the solution
            if (best_delta < 0)
            {
                int m = best_i;
                int n = best_k;
                while (m < n)
                {
                    int temp = current_solution[m];
                    current_solution[m] = current_solution[n];
                    current_solution[n] = temp;
                    m++;
                    n--;
                }
                improved = true;
            }

            // update the best solution if the current solution is better
            double current_distance = solution.getCost(current_solution, distances, numCities);
            double best_distance = solution.getCost(best_solution, distances, numCities);
            if (current_distance < best_distance)
            {
                for (int i = 0; i < numCities; i++)
                {
                    best_solution[i] = current_solution[i];
                }
            }
        }

        auto end_time = chrono::high_resolution_clock::now();
        // cout << "Time taken by function: " << chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << endl;
        // cout << "Number of iterations: " << iter_number << endl;
        // cout << "Number of solutions visited: " << solutions_visited << endl;
        return best_solution;
    }
};

int main()
{
    // read data from file
    string fileName = "tsp_instances/pr76.tsp"; // state file name
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
    int *best_solution;
    int time_limit = 2; // time limit in seconds
    time_t start_time = time(NULL);
    best_solution = algorithm.randomAlgorithm(sol, solution, distances, numCities, cost, start_time, time_limit);
    // print best solution
    for (int i = 0; i < numCities; i++)
    {
        cout << best_solution[i] << " ";
    }
    cout << endl;
    cout << "Cost: " << solution.getCost(best_solution, distances, numCities) << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // RANDOM WALK
    cout << endl
         << "RANDOM WALK" << endl;
    time_limit = 3; // time limit in seconds
    start_time = time(NULL);
    best_solution = algorithm.randomWalkAlgorithm(sol, solution, distances, numCities, cost, start_time, time_limit);
    // print best solution
    for (int i = 0; i < numCities; i++)
    {
        cout << best_solution[i] << " ";
    }
    cout << endl;
    cout << "Cost: " << solution.getCost(best_solution, distances, numCities) << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // NEAREST NEIGHBOR
    cout << endl
         << "NEAREST NEIGHBOR" << endl;
    // array for best solutions
    int **nn_solutions = new int *[10];
    for (int i = 0; i < 10; i++)
    {
        nn_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *nn_costs = new double[10];
    // perform nearest neighbor algorithm 10 times
    for (int i = 0; i < 10; i++)
    {
        nn_solutions[i] = algorithm.nearestNeighborAlgorithm(distances, numCities);
        nn_costs[i] = solution.getCost(nn_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    double avg_cost = 0;
    double min_cost = nn_costs[0];
    double max_cost = nn_costs[0];
    for (int i = 0; i < 10; i++)
    {
        avg_cost += nn_costs[i];
        if (nn_costs[i] < min_cost)
        {
            min_cost = nn_costs[i];
        }
        if (nn_costs[i] > max_cost)
        {
            max_cost = nn_costs[i];
        }
    }
    // save to file
    solution.saveToFile("nn_solutions.txt", nn_solutions, nn_costs, numCities);
    
    cout << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // GREEDY
    cout << endl
         << "GREEDY" << endl;
    // array for best solutions
    int **greedy_solutions = new int *[10];
    for (int i = 0; i < 10; i++)
    {
        greedy_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *greedy_costs = new double[10];
    // perform greedy algorithm 10 times
    for (int i = 0; i < 10; i++)
    {
        sol = solution.getRandomSolution(sol, numCities);
        greedy_solutions[i] = algorithm.greedyAlgorithm(sol, distances, numCities);
        greedy_costs[i] = solution.getCost(greedy_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    avg_cost = 0;
    min_cost = greedy_costs[0];
    max_cost = greedy_costs[0];
    for (int i = 0; i < 10; i++)
    {
        avg_cost += greedy_costs[i];
        if (greedy_costs[i] < min_cost)
        {
            min_cost = greedy_costs[i];
        }
        if (greedy_costs[i] > max_cost)
        {
            max_cost = greedy_costs[i];
        }
    }
    // save to file
    solution.saveToFile("greedy_solutions.txt", greedy_solutions, greedy_costs, numCities);
    cout << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // STEEPEST
    cout << endl
         << "STEEPEST" << endl;
    // array for best solutions
    int **steepest_solutions = new int *[10];
    for (int i = 0; i < 10; i++)
    {
        steepest_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *steepest_costs = new double[10];

    // perform steepest algorithm 10 times
    for (int i = 0; i < 10; i++)
    {
        sol = solution.getRandomSolution(sol, numCities);
        steepest_solutions[i] = algorithm.steepestAlgorithm(sol, distances, numCities);
        steepest_costs[i] = solution.getCost(steepest_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    avg_cost = 0;
    min_cost = steepest_costs[0];
    max_cost = steepest_costs[0];
    for (int i = 0; i < 10; i++)
    {
        avg_cost += steepest_costs[i];
        if (steepest_costs[i] < min_cost)
        {
            min_cost = steepest_costs[i];
        }
        if (steepest_costs[i] > max_cost)
        {
            max_cost = steepest_costs[i];
        }
    }
    // save to file
    solution.saveToFile("steepest_solutions.txt", steepest_solutions, steepest_costs, numCities);
    cout << endl;
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
