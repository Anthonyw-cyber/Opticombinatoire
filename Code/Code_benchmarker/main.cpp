


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>




using namespace std;

class Agent
{
public:
    vector<double> x;
    double B; // attractiveness of firefly
    double I; // intensity of firefly
    double objectivFunc;
    Agent();
    ~Agent();
};

double sphereFunction(const vector<double> &x);
double rastriginFunction(const vector<double> &x);
double ackleyFunction(const vector<double> &x);
double rosenbrockFunction(const vector<double> &x);

void generate_popultation(int m, int n, Agent firefly[], double x_r, double x_l);

double uRandom(double a);

int main()
{
    std::vector<int> dimensions = {10, 30, 50};
    std::vector<int> popSizes = {30, 50, 70};

    for (int dim : dimensions)
    {
        for (int popSize : popSizes)
        {
            cout << "dim: " << dim << " popSize: " << popSize << endl;
            string fileName = "results_dim_" + to_string(dim) + "_pop_" + to_string(popSize) + ".txt";
            ofstream outFile(fileName);

            if (!outFile)
            {
                cerr << "Cannot open output file: " << fileName << endl;
                return 1;
            }

            for (int run = 0; run < 10; ++run)
            {
                int pop_size = popSize;         // number of fireflies
                int n_dim = dim;             // dimension of the problem
                double lb = -10.0;      // lower bound
                double ub = 10.0;       // upper bound
                Agent firefly[pop_size];      // fireflies
                double alpha = 0.2;           // alpha, parameter of algorithm

                // initialization population of fireflies
                generate_popultation(pop_size, n_dim, firefly, ub, lb);

                double r = 0.0;        // distance between firefly
                double beta_base = 2;
                double gamma = 0.001;        // Gamma parameter
                double u = uRandom(alpha); // generate random u with uniform distribution
                int t = 1;
                double print;
                int epoch = 5000;
                double fitness = firefly[0].objectivFunc;

                while (t < epoch)
                {
                    for (size_t i = 0; i < pop_size; i++)
                    {
                        for (size_t j = 0; j < pop_size; j++)
                        {
                            for (size_t k = 0; k < n_dim; k++)
                            {
                                r += (firefly[i].x[k] - firefly[j].x[k]) * (firefly[i].x[k] - firefly[j].x[k]);
                            }
                            r = sqrt(r);

                            if (firefly[j].I > firefly[i].I)
                            {
                                firefly[i].B = beta_base * exp(-gamma * pow(r, 2));

                                for (size_t k = 0; k < n_dim; k++)
                                {
                                    u = uRandom(alpha);
                                    firefly[i].x[k] = firefly[i].x[k] + firefly[i].B * (firefly[j].x[k] - firefly[i].x[k]) + u;
                                }
                            }
                        }

                        // checking boundaries
                        for (size_t l = 0; l < n_dim; l++)
                        {
                            if (firefly[i].x[l] < lb)
                            {
                                firefly[i].x[l] = lb;
                            }
                            if (firefly[i].x[l] > ub)
                            {
                                firefly[i].x[l] = ub;
                            }
                        }

                        // new solutions and intensity
                        firefly[i].objectivFunc = rosenbrockFunction(firefly[i].x); // Change to the objective function you want to use
                        firefly[i].I = 1.0 / firefly[i].objectivFunc;
                    }

                    // search for the fitness
                    int best_index = 0;

                    for (size_t i = 0; i < pop_size; i++)
                    {
                        if (fitness > firefly[i].objectivFunc)
                        {
                            fitness = firefly[i].objectivFunc;
                            best_index = i;
                        }
                    }
                    for (size_t i = 0; i < n_dim; i++)
                    {
                        u = uRandom(alpha);
                        firefly[best_index].x[i] = firefly[best_index].x[i] + u;
                    }

                    t++;
                }
                outFile << fitness << endl;

            }

            outFile.close();
        }
    }

    return 0;
}

double sphereFunction(const vector<double> &x)
{
    double sum = 0.0;
    for (double val : x)
    {
        sum += val * val;
    }
    return sum;
}

double rastriginFunction(const vector<double> &x)
{
    double sum = 0.0;
    for (double val : x)
    {
        sum += val * val - 10 * cos(2 * M_PI * val) + 10;
    }
    return sum;
}

double ackleyFunction(const vector<double> &x)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (double val : x)
    {
        sum1 += val * val;
        sum2 += cos(2 * M_PI * val);
    }
    return -20 * exp(-0.2 * sqrt(sum1 / x.size())) - exp(sum2 / x.size()) + 20 + M_E;
}

double rosenbrockFunction(const vector<double> &x)
{
    double sum = 0.0;
    for (size_t i = 0; i < x.size() - 1; i++)
    {
        sum += 100 * pow(x[i + 1] - pow(x[i], 2), 2) + pow(1 - x[i], 2);
    }
    return sum;
}

void generate_popultation(int pop_size, int n_dim, Agent firefly[], double ub, double lb)
{
    // initialization fireflies
    for (size_t i = 0; i < pop_size; i++)
    {
        firefly[i].B = 1.0;
        double temp;
        for (size_t j = 0; j < n_dim; j++)
        {
            temp = (rand() / (double)RAND_MAX) * (ub - lb) + lb;
            firefly[i].x.push_back(temp);
        }
        firefly[i].objectivFunc = rosenbrockFunction(firefly[i].x); // Change to the objective function you want to use
        firefly[i].I = 1 / firefly[i].objectivFunc;
    }
}

double uRandom(double a)
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double E = distribution(generator);
    return a * (E - 0.5);
}

Agent::Agent()
{
}

Agent::~Agent()
{
}
