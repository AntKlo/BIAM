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

void inverse_order_of_subarray(int *curr_solution, int i, int k)
{
    while (i < k)
    {
        int temp = curr_solution[i];
        curr_solution[i] = curr_solution[k];
        curr_solution[k] = temp;
        i++;
        k--;
    }
}

void copy_array(int *to_set_arr, int *to_cp_arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        to_set_arr[i] = to_cp_arr[i];
    }
}

void print_avg_min_max_cost(double *costs, int size)
{
    int min = numeric_limits<int>::max();
    int max = numeric_limits<int>::min();
    int sum = 0;
    for (int i = 0; i < size; i++)
    {
        if (costs[i] < min)
        {
            min = costs[i];
        }
        if (costs[i] > max)
        {
            max = costs[i];
        }
        sum += costs[i];
    }
    cout << "Average cost: " << sum / size << endl;
    cout << "Minimum cost: " << min << endl;
    cout << "Maximum cost: " << max << endl;
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
        return numCities;
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
        return type;
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
        return cities;
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

    // function to calculate euclidean distances between cities
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

    // function to calculate distances between cities based on type
    double **calculateDistance(double **cities, double **distances, int numCities, string type)
    {
        if (type == "EUC_2D")
        {
            distances = euclideanDistance(cities, distances, numCities);
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
    int *randomAlgorithm(int *best_solution, Solution solution, double **distances, int numCities, double min_cost, int time_limit)
    {
        int iter_number = 0;
        auto start_time = chrono::high_resolution_clock::now();
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
        } while (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start_time).count() < time_limit);

        // cout << "Number of iterations: " << iter_number << endl;

        return best_solution;
    }
    // random walk algorithm
    int *randomWalkAlgorithm(int *best_solution, Solution solution, double **distances, int numCities, double min_cost, int time_limit)
    {
        int iter_number = 0;
        auto start_time = chrono::high_resolution_clock::now();
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
        } while (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start_time).count() < time_limit);        

        // cout << "Number of iterations: " << iter_number << endl;

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
    void greedyAlgorithm(int *initial_solution, int *new_solution, double **distances, int numCities)
    {
        Solution solution;
        auto start_time = chrono::high_resolution_clock::now();
        copy_array(new_solution, initial_solution, numCities);

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
                    double delta = round(100.0*(distances[new_solution[i - 1]][new_solution[k]] + distances[new_solution[i]][new_solution[k + 1]] - distances[new_solution[i - 1]][new_solution[i]] - distances[new_solution[k]][new_solution[k + 1]]))/100.0;
                    if (delta < 0)
                    {
                        inverse_order_of_subarray(new_solution, i, k);
                        improved = true;
                        break;
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
    }


    void steepestAlgorithm(int *initial_solution, int *new_solution, double **distances, int numCities)
    {
        Solution solution;
        auto start_time = chrono::high_resolution_clock::now();
        copy_array(new_solution, initial_solution, numCities);

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
                    double delta = round(100.0*(distances[new_solution[i - 1]][new_solution[k]] + distances[new_solution[i]][new_solution[k + 1]] - distances[new_solution[i - 1]][new_solution[i]] - distances[new_solution[k]][new_solution[k + 1]]))/100.0;
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
                inverse_order_of_subarray(new_solution, best_i, best_k);
                improved = true;
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        // cout << "Time taken by function: " << chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << endl;
        // cout << "Number of iterations: " << iter_number << endl;
        // cout << "Number of solutions visited: " << solutions_visited << endl;
    }
};

int main()
{
    // read data from file
    string instance_name = "ch150"; // instance name
    string fileName = "tsp_instances/"+instance_name+".tsp"; // state file name
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
    int number_of_iterations = 10;
    int *new_solution = new int[numCities];
    int *tmp_solution = new int[numCities];
    // RANDOM
    cout << endl
         << "RANDOM" << endl;
    // array for best solutions
    int **random_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        random_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *random_costs = new double[number_of_iterations];
    int time_limit = 10000000; // time limit in nanoseconds
    // perform algorithm number_of_iterations times
    for (int i = 0; i < number_of_iterations; i++)
    {
        // auto start_time = chrono::high_resolution_clock::now();
        random_solutions[i] = algorithm.randomAlgorithm(sol, solution, distances, numCities, cost, time_limit);
        random_costs[i] = solution.getCost(random_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(random_costs, number_of_iterations);

    // save to file
    solution.saveToFile("solution_random_"+instance_name+".txt", random_solutions, random_costs, numCities);
    delete[] random_solutions;
    delete[] random_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // RANDOM WALK
    cout << endl
         << "RANDOM WALK" << endl;
    time_limit = 10000000; // time limit in nanoseconds
    // array for best solutions
    int **rw_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        rw_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *rw_costs = new double[number_of_iterations];
    // perform algorithm number_of_iterations times
    for (int i = 0; i < number_of_iterations; i++)
    {
        rw_solutions[i] = algorithm.randomWalkAlgorithm(sol, solution, distances, numCities, cost, time_limit);
        rw_costs[i] = solution.getCost(rw_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(rw_costs, number_of_iterations);
    // save to file
    solution.saveToFile("solution_random_walk_"+instance_name+".txt", rw_solutions, rw_costs, numCities);
    delete[] rw_solutions;
    delete[] rw_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // NEAREST NEIGHBOR
    cout << endl
         << "NEAREST NEIGHBOR" << endl;
    // array for best solutions
    int **nn_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        nn_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *nn_costs = new double[number_of_iterations];
    // perform nearest neighbor algorithm number_of_iterations times
    for (int i = 0; i < number_of_iterations; i++)
    {
        nn_solutions[i] = algorithm.nearestNeighborAlgorithm(distances, numCities);
        nn_costs[i] = solution.getCost(nn_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(nn_costs, number_of_iterations);
    // save to file
    solution.saveToFile("solution_nn_"+instance_name+".txt", nn_solutions, nn_costs, numCities);
    delete[] nn_solutions;
    delete[] nn_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // GREEDY
    cout << endl
         << "GREEDY" << endl;
    // array for best solutions
    int **greedy_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        greedy_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *greedy_costs = new double[number_of_iterations];
    // perform greedy algorithm number_of_iterations times
    for (int i = 0; i < number_of_iterations; i++)
    {
        sol = solution.getRandomSolution(sol, numCities);
        algorithm.greedyAlgorithm(sol, new_solution, distances, numCities);
        copy_array(greedy_solutions[i], new_solution, numCities);
        greedy_costs[i] = solution.getCost(greedy_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(greedy_costs, number_of_iterations);
    // save to file
    solution.saveToFile("solution_greedy_"+instance_name+".txt", greedy_solutions, greedy_costs, numCities);
    delete[] greedy_solutions;
    delete[] greedy_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------
    // STEEPEST
    cout << endl
         << "STEEPEST" << endl;
    // array for best solutions
    int **steepest_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        steepest_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *steepest_costs = new double[number_of_iterations];

    // perform steepest algorithm number_of_iterations times
    for (int i = 0; i < number_of_iterations; i++)
    {
        sol = solution.getRandomSolution(sol, numCities);
        algorithm.steepestAlgorithm(sol, new_solution, distances, numCities);
        steepest_solutions[i] = new_solution;
        steepest_costs[i] = solution.getCost(steepest_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(steepest_costs, number_of_iterations);
    // save to file
    solution.saveToFile("solution_steepest_"+instance_name+".txt", steepest_solutions, steepest_costs, numCities);
    delete[] steepest_solutions;
    delete[] steepest_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;

    // free the dynamically allocated memory
    for (int i = 0; i <= numCities; i++)
    {
        delete[] cities[i];
        // delete[] distances[i];  // in my case Exception has occurred Segmentation fault -> just comment for now to make it work.
    }
    delete[] cities;
    delete[] distances;
    delete[] sol;
    delete[] new_solution;
    delete[] tmp_solution;
    return 0;
}
