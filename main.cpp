#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <typeinfo>

using namespace std;

void swap(int a, int b)
{
    /**
     * Swap two numbers without using a third variable.
     *
     */
    b = a + b;
    a = b - a;
    b = b - a;
}

class Problem
{
public:
    // function to randomize vector
    void randomize(int arr[], int n)
    {
        // Use a different seed value so that
        // we don't get same result each time
        // we run this program
        // srand(time(NULL));

        // Start from the last element and swap
        // one by one. We don't need to run for
        // the first element that's why i > 0
        for (int i = n - 1; i > 0; i--)
        {
            // Pick a random index from 0 to i
            int j = rand() % (i + 1);

            // Swap arr[i] with the element
            // at random index
            swap(arr[i], arr[j]);
        }
    }

    // function that calculates euclidean distance between two points
    double euclideanDistance(double x1, double y1, double x2, double y2)
    {
        return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    }   


    double** readFromFile(string fileName)
    {
        int numCities = 0;
        double **cities;
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
            double **cities = new double *[numCities];  

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
                    string x = str_copy.substr(0);
                    pos = str_copy.find(" ");
                    str_copy = str_copy.substr(pos + 1);
                    string y = str_copy.substr(0);
                    cities[i][0] = stod(x);
                    cities[i][1] = stod(y);
                    cout << cities[i][0] << "   " << cities[i][1] << endl;
                }
            }
        }
        // print length of cities
        cout << "Length of cities: " << sizeof(cities) << endl;
        // print all the coordinates
        for (int i = 0; i < numCities; i++)
        {
            cout << cities[i][0] << " " << cities[i][1] << endl;
        }
        return cities;
    }
};

int main()
{
    Problem problem;
    double ** cities;
    cities = problem.readFromFile("tsp_instances/ulysses22.tsp");
    return 0;
}