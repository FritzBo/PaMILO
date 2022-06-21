/**
 * @file dual_benson_scalarizer.h
 * @author Fritz Bökler
 * @brief
 * @date 27.09.2013
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

//#define WEPS
//#define VALEPS
//#define DMEAS

#pragma once

#include <functional>
#include <vector>

#include <pamilo/basic/point.h>

// todo: delete this include
// #include <iomanip>

namespace pamilo {

/**
 * @brief Objects of this class run the Dual Benson Algorithm
 *
 * @tparam OnlineVertexEnumerator Vertex enumerator for the vertex enumeration. Any class deriving
 * from AbstractOnlineVertexEnumerator is suited.
 * @tparam SolType Type to store the solution in. String by default.
 */
template <typename OnlineVertexEnumerator, typename SolType = std::string>
class DualBensonScalarizer
{
public:
    /**
     * @brief Construct a new Dual Benson Scalarizer and initializes its parameters.
     *
     * @param solver Callable object that solves weighted sum problems (for example
     * ILPSolverAdaptor)
     * @param printSol Callable object which stores solutions in the form of pair<SolType, Point*>
     * (for example ILPSolverPrinter)
     * @param dimension dimension of objective space
     * @param epsilon epsilon value for floating point comparisons
     * @param pepsilon epsilon value for point equality comparisons
     * @param veEpsilon epsilon value for floating points comparisons in vertex enumeration
     */
    DualBensonScalarizer(
        std::function<double(const Point &weighting, Point &value, SolType &sol)> solver,
        std::function<void(std::pair<SolType, Point *>, bool)> printSol, unsigned int dimension,
        double epsilon, double pEpsilon, double veEpsilon)
        : dimension_(dimension)
        , epsilon_(epsilon)
        , pEpsilon_(pEpsilon)
        , veEpsilon_(veEpsilon)
        , solver_(solver)
        , printSol_(printSol)
        , vertex_container(nullptr)
        , vertices_(0)
        , facets_(0)
        , ve_time(0)
    {
        if (veEpsilon_ == -1)
        {
            veEpsilon_ = epsilon_;
        }
    }

    /**
     * @brief Destructor
     *
     */
    ~DualBensonScalarizer()
    {
        if (vertex_container)
        {
            delete vertex_container;
        }
    }

    /**
     * @brief Calculates the non dominated extreme points and stores them in solutions.
     * Additionally, every solutions is input into the callable object printSol given to the
     * constructor.
     *
     * @param solutions array of pairs onto which all solutions are added
     */
    void Calculate_solutions(std::vector<std::pair<SolType, Point *>> &solutions);

    /**
     * @brief Returns cpu time used for vertex enumeration since instantiation of this object
     *
     * @return double
     */
    double vertex_enumeration_time();

    /**
     * @brief number of vertices of the upper image found. Initially 0
     *
     * @return int
     */
    int number_vertices()
    {
        return vertices_;
    }

    /**
     * @brief number of facets of the upper image found. Initially 0
     *
     * @return int
     */
    int number_facets()
    {
        return facets_;
    }

protected:
    /**
     * @brief number of objective space dimensions
     *
     */
    unsigned int dimension_;

    /**
     * @brief epsilon value for floating point comparisons
     *
     */
    double epsilon_;

    /**
     * @brief epsilon value for euclidean distances of points comparisons
     *
     */
    double pEpsilon_;

    /**
     * @brief epsilon value for floating points comparisons in vertex enumeration
     *
     */
    double veEpsilon_;

    /**
     * @brief Callable object that solves weighted sum problems (for example ILPSolverAdaptor)
     *
     */
    std::function<double(const Point &, Point &, SolType &sol)> solver_;

    /**
     * @brief
     *
     */
    std::function<void(std::pair<SolType, Point *>, bool)> printSol_;

private:
    /**
     * @brief vertex enumerator object
     *
     */
    OnlineVertexEnumerator *vertex_container;

    /**
     * @brief number of vetices of the upper image
     *
     */
    int vertices_;

    /**
     * @brief number of facets of the upper image
     *
     */
    int facets_;

    /**
     * @brief cpu time for vertex enumeration
     *
     */
    double ve_time;
};

template <typename OnlineVertexEnumerator, typename SolType>
void DualBensonScalarizer<OnlineVertexEnumerator, SolType>::Calculate_solutions(
    std::vector<std::pair<SolType, Point *>> &solutions)
{
    vertices_ = 1;
    int iteration_counter = 0;
    facets_ = 1;

    Point v(dimension_);
    Point value(dimension_);

    for (unsigned int i = 0; i < dimension_ - 1; ++i)
        v[i] = 0;
    v[0] = 1;

    SolType sol;
    solver_(v, value, sol);
    auto solPair = std::make_pair(sol, new Point(value));
    solutions.push_back(solPair);
    printSol_(solPair, true);

    clock_t start = clock();
    if (vertex_container)
    {
        delete vertex_container;
    }
    vertex_container = new OnlineVertexEnumerator(value, dimension_, veEpsilon_);
    delete vertex_container->next_vertex();
    ve_time += (clock() - start) / (double)CLOCKS_PER_SEC;

    Point *candidate;
    Point weighting(dimension_);
    Point inequality(dimension_);

    double scalar_value;
    while (vertex_container->has_next())
    {
        iteration_counter++;

#ifndef NDEBUG
        std::cout << "Iteration: " << iteration_counter << std::endl;
#endif

        candidate = vertex_container->next_vertex();

#ifndef NDEBUG
        std::cout << "New candidate: " << *candidate << std::endl;
#endif

        double sum = 0;
        for (unsigned int i = 0; i < dimension_ - 1; ++i)
        {
            weighting[i] = (*candidate)[i];
            sum += (*candidate)[i];
        }
        weighting[dimension_ - 1] = 1 - sum;

        for (unsigned int i = 0; i < dimension_; ++i)
        {
            value[i] = 0;
            if (weighting[i] < 0)
            {
#ifndef NDEBUG
                std::printf("\n\n\nweight eps: %E\n\n\n", weighting[i]);
#endif
#ifdef WEPS
                weighting[i] = 0;
#endif
            }
        }

#ifndef NDEBUG
        std::cout << "weighting: " << weighting << std::endl;
#endif

        // todo delete
        //  std::cout << std::setprecision(std::numeric_limits<double>::digits10) << "\nWeighting: "
        //  << weighting << "\tScalar: " << candidate->operator[](dimension_-1) << "\n";

        SolType sol;
        scalar_value = solver_(weighting, value, sol);

        bool nextToExisting = false;
        if (pEpsilon_ >= 0)
        {
            nextToExisting = (std::find_if(solutions.begin(), solutions.end(),
                                           [value, this](const std::pair<SolType, Point *> &s) {
                                               Point dif = *(s.second) - value;
                                               return sqrt(dif * dif) <= this->pEpsilon_;
                                           }) != solutions.end());
        }

#ifndef NDEBUG
        std::cout << "scalar value: " << scalar_value << std::endl;
        std::cout << "value vector: " << value << std::endl;
        std::cout << "candidate value: " << (*candidate)[dimension_ - 1] << std::endl;
#endif
        for (unsigned int i = 0; i < dimension_ - 1; ++i)
        {
            double val = value[i] - value[dimension_ - 1];
            if (abs(val) < epsilon_)
            {
#ifndef NDEBUG
                std::cout << "\n\n\nval eps. val: " << val << ", eps: " << epsilon_
                          << " component: " << i << ", from value[i] = " << value[i]
                          << " and value[dim-1] = " << value[dimension_ - 1] << "\n\n\n";
#endif
#ifdef VALEPS
                val = 0;
#endif
            }
            inequality[i] = val;
        }
        inequality[dimension_ - 1] = -1;

        if (scalar_value - (*candidate)[dimension_ - 1] > -epsilon_ || nextToExisting
#ifdef DMEAS
            || vertex_container->getDistance(*candidate, inequality, -value[dimension_ - 1]) >
                   -epsilon_
#endif
        )
        {
            // todo delete
            //  std::cout << std::setprecision(std::numeric_limits<double>::digits10) << "Facet
            //  verified\tScalar " << scalar_value << "\tDif: " << scalar_value -
            //  (*candidate)[dimension_ - 1]  << std::endl;
            facets_++;
#ifndef NDEBUG
            std::cout << "found a new permanent extreme point. continuing." << std::endl;
#endif
        }
        else
        {
            // todo delete
            //  std::cout << std::setprecision(std::numeric_limits<double>::digits10) <<  "Point: "
            //  << value << "\tScalar: " << scalar_value << "\tDif: " << scalar_value -
            //  (*candidate)[dimension_ - 1] << std::endl;

            // todo delete
            // if (nextToExisting) {
            //     std::cout << "Prune" << std::endl;
            // }
            // std::cout << value << std::endl;

            clock_t start = clock();
            vertex_container->add_hyperplane(*candidate, inequality, -value[dimension_ - 1]);
            ve_time += (clock() - start) / (double)CLOCKS_PER_SEC;
            vertices_++;

            auto solPair = std::make_pair(sol, new Point(value));
            solutions.push_back(solPair);
            printSol_(solPair, false);
        }

        delete candidate;
    }

#ifndef NDEBUG
    std::cout << "Found " << facets_ << " nondominated value vectors in " << iteration_counter
              << " iterations." << std::endl;
    std::cout << "Where " << vertices_ << " weightings have been explored." << std::endl;
#else
    std::cout << "Found " << vertices_ << " non-dominated extreme points." << std::endl;
#endif
}

template <typename OnlineVertexEnumerator, typename SolType>
double DualBensonScalarizer<OnlineVertexEnumerator, SolType>::vertex_enumeration_time()
{
    return ve_time;
}
}  // namespace pamilo
