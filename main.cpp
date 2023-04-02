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
#include <random>

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

int uniform_random(int min, int max)
{
    static bool first = true;
    if (first)
    {
        srand(time(NULL));
        first = false;
    }
    int N = max - min + 1;
    return min + rand() / (RAND_MAX / N + 1);
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
    int *getInitialSolution(int numCities)
    {
        int *sol = new int[numCities + 1];
        for (int i = 0; i < numCities; i++)
        {
            sol[i] = i;
        }
        sol[numCities] = sol[0]; // add the first city to the end of the array to allow easier computations
        return sol;
    }
    // function to generate a random solution
    void makeRandom(int *arr, int size)
    {
        // Iterate through the array starting from the last element
        for (int i = size - 1; i > 0; i--)
        {
            // Generate a random index between 0 and i (inclusive)
            int j = uniform_random(0, i);
            // Swap the current element with a random element before it
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
        arr[size] = arr[0];
    }

    void check_solution(int *sol, int size)
    {
        // checks whether there exist pair of cities that are repeated in the solution
        for (int i = 0; i < size; i++)
        {
            for (int j = i + 1; j < size; j++)
            {
                if (sol[i] == sol[j])
                {
                    cout << "ERROR: repeated city in solution" << endl;
                    exit(1);
                }
            }
        }
        if (sol[size] != sol[0])
        {
            cout << "ERROR: ending city not the same as the beginning one" << endl;
            exit(1);
        }
        std::cout << "Checked! All fine!" << std::endl;
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
        for (int i = 0; i < numCities; i++) // array should have first city at the end as well
        {
            cost += distances[sol[i]][sol[i + 1]];
        }
        return cost;
    }
    double getCostWithoutAppended(int *sol, double **distances, int numCities)
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
    void randomAlgorithm(int *best_solution, int *current_solution, Solution solution_utilities, double **distances, int numCities, double min_cost, int time_limit)
    {
        int iter_number = 0;
        auto start_time = chrono::high_resolution_clock::now();
        double new_cost;
        do
        {
            iter_number++;
            solution_utilities.makeRandom(current_solution, numCities); // isn't it the same as initial solution?
            new_cost = solution_utilities.getCost(current_solution, distances, numCities);
            if (new_cost < min_cost)
            {
                min_cost = new_cost;
                copy_array(best_solution, current_solution, numCities);
            }
        } while (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start_time).count() < time_limit);

        cout << "Number of iterations: " << iter_number << endl;
    }

    void randomWalkAlgorithm(int *best_solution, int *current_solution, Solution solution_utilities, double **distances, int numCities, int time_limit, int **edge_pairs_nonadjacent, int edge_pairs_nonadjacent_SIZE)
    {
        int iter_number = 0;
        auto start_time = chrono::high_resolution_clock::now();
        int edges_pair_id;
        double delta;
        int q, w;
        double delta_sum = 0.0;
        copy_array(best_solution, current_solution, numCities); // in case when there is no improvement first solution should be set as the best
        do
        {
            iter_number++;
            edges_pair_id = uniform_random(0, edge_pairs_nonadjacent_SIZE - 1);
            q = edge_pairs_nonadjacent[edges_pair_id][0];
            w = edge_pairs_nonadjacent[edges_pair_id][1];
            delta = distances[current_solution[q - 1]][current_solution[w]] + distances[current_solution[q]][current_solution[w + 1]] - distances[current_solution[q - 1]][current_solution[q]] - distances[current_solution[w]][current_solution[w + 1]];
            delta_sum += delta;
            inverse_order_of_subarray(current_solution, q, w);
            if (round(delta_sum * 1000) < 0)
            {
                copy_array(best_solution, current_solution, numCities + 1);
                delta_sum = 0.0;
            }
        } while (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start_time).count() < time_limit);
        cout << "Number of iterations: " << iter_number << endl;
    }

    // nearest neighbor algorithm
    int *nearestNeighborAlgorithm(double **distances, int numCities)
    {
        auto start_time = chrono::high_resolution_clock::now();
        int *new_solution = new int[numCities];
        new_solution[0] = uniform_random(0, numCities - 1);
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

    void greedyAlgorithm(int *new_solution, double **distances, int numCities)
    {
        Solution solution;
        bool improved = true;
        int iter_number = 0;
        int solutions_visited = 0;
        int i;
        double delta;
        while (improved)
        {
            iter_number++;
            improved = false;
            for (int i = 1; i < numCities - 1; i++)
            {
                if (i == 1)
                { // we don't use edge between 0 and numCities - 1 because it is shares node 0
                    for (int k = i + 1; k < numCities - 1; k++)
                    {
                        solutions_visited++;
                        delta = distances[new_solution[i - 1]][new_solution[k]] + distances[new_solution[i]][new_solution[k + 1]] - distances[new_solution[i - 1]][new_solution[i]] - distances[new_solution[k]][new_solution[k + 1]];
                        if (round(delta * 1000) < 0)
                        {
                            inverse_order_of_subarray(new_solution, i, k);
                            improved = true;
                            break;
                        }
                    }
                }
                else
                {
                    for (int k = i + 1; k < numCities; k++) // here we use edge between 0 and numCities - 1, which was added at the end of the solution
                    {
                        solutions_visited++;
                        delta = distances[new_solution[i - 1]][new_solution[k]] + distances[new_solution[i]][new_solution[k + 1]] - distances[new_solution[i - 1]][new_solution[i]] - distances[new_solution[k]][new_solution[k + 1]];
                        if (round(delta * 1000) < 0)
                        {
                            inverse_order_of_subarray(new_solution, i, k);
                            improved = true;
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
    }

    void steepestAlgorithm(int *new_solution, double **distances, int numCities)
    {
        Solution solution;
        auto start_time = chrono::high_resolution_clock::now();
        bool improved = true;
        int iter_number = 0;
        int solutions_visited = 0;
        double best_delta;
        int best_i;
        int best_k;
        double delta;
        while (improved)
        {
            improved = false;
            iter_number++;
            best_delta = 0.0;
            best_i = -1;
            best_k = -1;
            // iterate over all possible moves
            for (int i = 1; i < numCities - 1; i++)
            {
                if (i == 1)
                { // we don't use edge between 0 and numCities - 1 because it is shares node 0
                    for (int k = i + 1; k < numCities - 1; k++)
                    {
                        solutions_visited++;
                        delta = distances[new_solution[i - 1]][new_solution[k]] + distances[new_solution[i]][new_solution[k + 1]] - distances[new_solution[i - 1]][new_solution[i]] - distances[new_solution[k]][new_solution[k + 1]];
                        if (delta < best_delta)
                        {
                            best_delta = delta;
                            best_i = i;
                            best_k = k;
                        }
                    }
                }
                else
                {
                    for (int k = i + 1; k < numCities; k++) // here we use edge between 0 and numCities - 1, which was added at the end of the solution
                    {
                        solutions_visited++;
                        delta = distances[new_solution[i - 1]][new_solution[k]] + distances[new_solution[i]][new_solution[k + 1]] - distances[new_solution[i - 1]][new_solution[i]] - distances[new_solution[k]][new_solution[k + 1]];
                        if (delta < best_delta)
                        {
                            best_delta = delta;
                            best_i = i;
                            best_k = k;
                        }
                    }
                }
            }

            // apply the best move if it improves the solution
            if (best_delta < 0.0)
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

    void simulatedAnnealingAlgorithm(int *new_solution, double **distances, int numCities, double initialTemperature, double endTemperature, int length, double alpha, int **edge_pairs_nonadjacent, int edge_pairs_nonadjacent_SIZE)
    {
        double temperature = initialTemperature;
        double delta;
        int edges_pair_id;
        int q;
        int w;
        int non_improving_iterations = 0;
        while (temperature > endTemperature && non_improving_iterations < 10 * length)
        {
            for (int l = 0; l < length; l++)
            {
                // get neighbour solution
                edges_pair_id = uniform_random(0, edge_pairs_nonadjacent_SIZE - 1);
                q = edge_pairs_nonadjacent[edges_pair_id][0];
                w = edge_pairs_nonadjacent[edges_pair_id][1];
                // calculate delta
                delta = distances[new_solution[q - 1]][new_solution[w]] + distances[new_solution[q]][new_solution[w + 1]] - distances[new_solution[q - 1]][new_solution[q]] - distances[new_solution[w]][new_solution[w + 1]];
                // accept or reject neighbour solution
                if (round(delta * 1000) < 0)
                {
                    // accept neighbour solution
                    inverse_order_of_subarray(new_solution, q, w);
                    non_improving_iterations = 0;
                }
                else
                {
                    double probability = exp(-delta / temperature);
                    double random = (double)rand() / RAND_MAX;
                    if (random < probability)
                    {
                        // accept neighbour solution
                        inverse_order_of_subarray(new_solution, q, w);
                        non_improving_iterations = 0;
                    }
                    else
                    {
                        // reject neighbour solution
                        non_improving_iterations++;
                    }
                }
            }
            // decrease temperature
            temperature = temperature * alpha;
        }
    }
};

void experiment_random(int *sol, Solution solution_utilities, int numCities, double time_limit, double **distances, Algorithm algorithm_functions, string instance_name)
{
    int number_of_iterations = 10;
    cout << endl
         << "RANDOM" << endl;
    // array for best solutions
    int **random_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        random_solutions[i] = new int[numCities + 1];
    }
    // array for costs of best solutions
    double *random_costs = new double[10];
    // perform algorithm 10 times
    int cost;
    for (int i = 0; i < 10; i++)
    {
        solution_utilities.makeRandom(sol, numCities);
        cost = solution_utilities.getCost(sol, distances, numCities);
        int *best_solution = new int[numCities + 1];
        algorithm_functions.randomAlgorithm(best_solution, sol, solution_utilities, distances, numCities, cost, time_limit);
        random_solutions[i] = best_solution;
        random_costs[i] = solution_utilities.getCost(random_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(random_costs, number_of_iterations);

    // save to file
    solution_utilities.saveToFile("solution_random_" + instance_name + ".txt", random_solutions, random_costs, numCities);
    delete[] random_solutions;
    delete[] random_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;
}

int **make_randomWalk_array_of_edge_pairs_non_adjacent(int numCities, Solution solution_utilities)
{
    int *sol = solution_utilities.getInitialSolution(numCities);
    int **edge_pairs = new int *[(numCities - 3) * numCities / 2]; // equivalent to: [(numCities - 2) * (numCities - 3) / 2 + (numCities - 3)];
    int counter = 0;
    for (int i = 1; i < numCities - 1; i++)
    {
        if (i == 1)
        {
            for (int j = i + 1; j < numCities - 1; j++)
            {
                edge_pairs[counter] = new int[2];
                edge_pairs[counter][0] = sol[i];
                edge_pairs[counter][1] = sol[j];
                counter++;
            }
        }
        else
        {
            for (int j = i + 1; j < numCities; j++)
            {
                edge_pairs[counter] = new int[2];
                edge_pairs[counter][0] = sol[i];
                edge_pairs[counter][1] = sol[j];
                counter++;
            }
        }
    }
    delete[] sol;
    return edge_pairs;
}

void remove_randomWalk_array_of_edges(int **edge_pairs, int numCities)
{
    for (int i = 0; i < (numCities - 3) * numCities / 2; i++)
    {
        delete[] edge_pairs[i];
    }
    delete[] edge_pairs;
}

void experiment_randomWalk(int *sol, Solution solution_utilities, int numCities, double time_limit, double **distances, Algorithm algorithm_functions, string instance_name)
{
    int number_of_iterations = 10;
    cout << endl
         << "RANDOM WALK" << endl;
    solution_utilities.makeRandom(sol, numCities);
    int **edge_pairs = make_randomWalk_array_of_edge_pairs_non_adjacent(numCities, solution_utilities);
    int edge_pairs_nonadjacent_size = (numCities - 3) * numCities / 2;
    // array for best solutions
    int **rw_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        rw_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *rw_costs = new double[10];
    int cost;
    int *new_solution = new int[numCities];
    int numIterations = 10;
    for (int i = 0; i < numIterations; i++)
    {
        solution_utilities.makeRandom(sol, numCities);
        algorithm_functions.randomWalkAlgorithm(new_solution, sol, solution_utilities, distances, numCities, time_limit, edge_pairs, edge_pairs_nonadjacent_size);
        rw_solutions[i] = new_solution;
        rw_costs[i] = solution_utilities.getCost(rw_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(rw_costs, number_of_iterations);
    // save to file
    solution_utilities.saveToFile("solution_random_walk_" + instance_name + ".txt", rw_solutions, rw_costs, numCities);
    delete[] rw_solutions;
    delete[] rw_costs;
    delete[] new_solution;
    cout << endl;
    cout << "-----------------------------------------" << endl;
    remove_randomWalk_array_of_edges(edge_pairs, numCities);
}

void experiment_nearest_neighbor(Solution solution_utilities, int numCities, double **distances, Algorithm algorithm_functions, string instance_name)
{
    int number_of_iterations = 10;
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
        nn_solutions[i] = algorithm_functions.nearestNeighborAlgorithm(distances, numCities);
        nn_costs[i] = solution_utilities.getCostWithoutAppended(nn_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(nn_costs, number_of_iterations);
    // save to file
    solution_utilities.saveToFile("solution_nn_" + instance_name + ".txt", nn_solutions, nn_costs, numCities);
    delete[] nn_solutions;
    delete[] nn_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;
}

void experiment_greedy(int *sol, Solution solution_utilities, int numCities, double **distances, Algorithm algorithm_functions, string instance_name)
{
    int number_of_iterations = 10;
    cout << endl
         << "GREEDY" << endl;
    // array for best solutions
    int **greedy_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        greedy_solutions[i] = new int[numCities + 1];
    }
    // array for costs of best solutions
    double *greedy_costs = new double[10];
    // perform greedy algorithm 10 times
    for (int i = 0; i < number_of_iterations; i++)
    {
        solution_utilities.makeRandom(sol, numCities);
        algorithm_functions.greedyAlgorithm(sol, distances, numCities); // works in situ on 'sol' variable
        copy_array(greedy_solutions[i], sol, numCities + 1);
        greedy_costs[i] = solution_utilities.getCost(greedy_solutions[i], distances, numCities);
    }
    print_avg_min_max_cost(greedy_costs, number_of_iterations);
    solution_utilities.saveToFile("solution_greedy_" + instance_name + ".txt", greedy_solutions, greedy_costs, numCities);
    delete[] greedy_solutions;
    delete[] greedy_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;
}

void experiment_steepest(int *sol, Solution solution_utilities, int numCities, double **distances, Algorithm algorithm_functions, string instance_name)
{
    int number_of_iterations = 10;
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
        solution_utilities.makeRandom(sol, numCities);
        algorithm_functions.steepestAlgorithm(sol, distances, numCities); // works in situ on 'sol' variable
        copy_array(steepest_solutions[i], sol, numCities + 1);
        steepest_costs[i] = solution_utilities.getCost(steepest_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(steepest_costs, number_of_iterations);
    // save to file
    solution_utilities.saveToFile("solution_steepest_" + instance_name + ".txt", steepest_solutions, steepest_costs, numCities);
    delete[] steepest_solutions;
    delete[] steepest_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;
}

void experiment_simulated_annealing(int *sol, Solution solution_utilities, int numCities, double **distances, Algorithm algorithm_functions, string instance_name, double initial_temperature, double endTemperature, int length, double alpha)
{
    int number_of_iterations = 10;
    cout << endl
         << "SIMULATED ANNEALING" << endl;
    int **edge_pairs = make_randomWalk_array_of_edge_pairs_non_adjacent(numCities, solution_utilities);
    int edge_pairs_nonadjacent_size = (numCities - 3) * numCities / 2;
    // array for best solutions
    int **sa_solutions = new int *[number_of_iterations];
    for (int i = 0; i < number_of_iterations; i++)
    {
        sa_solutions[i] = new int[numCities];
    }
    // array for costs of best solutions
    double *sa_costs = new double[number_of_iterations];

    // perform simulated annealing algorithm number_of_iterations times
    for (int i = 0; i < number_of_iterations; i++)
    {
        solution_utilities.makeRandom(sol, numCities);
        algorithm_functions.simulatedAnnealingAlgorithm(sol, distances, numCities, initial_temperature, endTemperature, length, alpha, edge_pairs, edge_pairs_nonadjacent_size); // works in situ on 'sol' variable
        copy_array(sa_solutions[i], sol, numCities + 1);
        sa_costs[i] = solution_utilities.getCost(sa_solutions[i], distances, numCities);
    }
    // print average, min and max cost
    print_avg_min_max_cost(sa_costs, number_of_iterations);
    // save to file
    solution_utilities.saveToFile("solution_sa_" + instance_name + ".txt", sa_solutions, sa_costs, numCities);
    delete[] sa_solutions;
    delete[] sa_costs;
    cout << endl;
    cout << "-----------------------------------------" << endl;
}

int main()
{
    srand(time(NULL));
    // read data from file
    string fileName = "tsp_instances/pr76.tsp"; // state file name
    string instance_name = "pr76";
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
    Solution solution_utilities;
    distances = solution_utilities.calculateDistance(cities, distances, numCities, type); // calculate distance between cities

    int *sol = solution_utilities.getInitialSolution(numCities);
    solution_utilities.makeRandom(sol, numCities);                       // shuffles solution to produce random solution
    double cost = solution_utilities.getCost(sol, distances, numCities); // get cost of initial solution
    cout << "Initial solution cost: " << cost << endl;
    cout << "-----------------------------------------" << endl;

    // -----------------------------------------

    // perform different algorithms
    Algorithm algorithm;
    // int time_limit = 3000000; // time limit in nanoseconds; used for random, randomWalk

    // experiment_random(sol, solution_utilities, numCities, time_limit, distances, algorithm, instance_name);
    // experiment_randomWalk(sol, solution_utilities, numCities, time_limit, distances, algorithm, instance_name);
    // experiment_nearest_neighbor(solution_utilities, numCities, distances, algorithm, instance_name);
    // experiment_greedy(sol, solution_utilities, numCities, distances, algorithm, instance_name);
    // experiment_steepest(sol, solution_utilities, numCities, distances, algorithm, instance_name);
    experiment_simulated_annealing(sol, solution_utilities, numCities, distances, algorithm, instance_name, 1000, 0.01, 100, 0.99);

    // free the dynamically allocated memory
    for (int i = 0; i <= numCities; i++)
    {
        delete[] cities[i];
        // delete[] distances[i];  // in my case Exception has occurred Segmentation fault -> just comment for now to make it work.
    }
    delete[] cities;
    delete[] distances;
    delete[] sol;
    return 0;
}
