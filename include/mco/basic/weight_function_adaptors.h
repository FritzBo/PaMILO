//
//  LexWeightFunctionAdaptor.h
//  mco
//
//  Created by Fritz Bökler on 04.04.14.
//
//

#ifndef mco_LexWeightFunctionAdaptor_h
#define mco_LexWeightFunctionAdaptor_h

#include <functional>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {
    
class WeightFunctionAdaptor {
public:
    WeightFunctionAdaptor(std::function<const Point *(const ogdf::edge)> cost_function,
                          const Point& weighting)
    :   weighting_(weighting),
    cost_function_(cost_function) { }
    
    double operator()(ogdf::edge e) {
        return *cost_function_(e) * weighting_;
    }
    
private:
    const Point& weighting_;
    std::function<const Point*(const ogdf::edge)> cost_function_;
    
};
    
    
    
    
class LexWeightFunctionAdaptor {
public:
    LexWeightFunctionAdaptor(const ogdf::Graph& graph,
                             std::function<const Point* (const ogdf::edge)> cost_function,
                             const Point& weighting)
    :   weighting_(weighting),
    cost_function_(cost_function),
    weighted_costs_(graph, nullptr),
    target_dimension_(weighting.dimension() + 1) { }
    
    ~LexWeightFunctionAdaptor();
    
    inline Point* operator()(ogdf::edge e);
    
private:
    const Point& weighting_;
    std::function<const Point*(const ogdf::edge)> cost_function_;
    ogdf::EdgeArray<Point *> weighted_costs_;
    unsigned target_dimension_;
    
};
    
inline LexWeightFunctionAdaptor::
~LexWeightFunctionAdaptor() {
    for(auto point : weighted_costs_) {
        delete point;
    }
}
    
inline Point* LexWeightFunctionAdaptor::
operator()(const ogdf::edge e) {
    
    if(weighted_costs_[e] != nullptr) {
        return weighted_costs_(e);
    }
    
    auto new_point = new Point(target_dimension_);
    const Point& cost = *cost_function_(e);

    std::copy(cost.cbegin(), cost.cend(), new_point->begin() + 1);
    new_point->operator[](0) = weighting_ * cost;
    
    weighted_costs_(e) = new_point;
    
    return new_point;
}

    
}

#endif
