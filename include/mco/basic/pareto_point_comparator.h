//
//  pareto_point_comparator.h
//  mco
//
//  Created by Fritz Bökler on 31.03.14.
//
//

#ifndef mco_pareto_point_comparator_h
#define mco_pareto_point_comparator_h

class ParetoDominationPointComparator {
public:
    inline ParetoDominationPointComparator(double epsilon = 0)
    : epsilon_(epsilon) { }
    
    inline bool operator()(const Point * p1,
                           const Point * p2) const noexcept;
    
    inline bool operator()(const Point& p1,
                           const Point& p2) const noexcept;
    
    inline static bool dominates(const Point& p1,
                                 const Point& p2,
                                 double epsilon) noexcept;
    
private:
    const double epsilon_;
};

inline bool ParetoDominationPointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {
    
    return dominates(p1, p2, epsilon_);
}

inline bool ParetoDominationPointComparator::
operator()(const Point * p1,
           const Point * p2) const noexcept {
    
    return dominates(*p1, *p2, epsilon_);
}

inline bool ParetoDominationPointComparator::
dominates(const Point &p1,
          const Point &p2,
          double epsilon) noexcept {
    
    assert(p1.dimension() == p2.dimension());
    
    unsigned dimension = p1.dimension();
    
    bool equal = true;
    
    for(unsigned int i = 0; i < dimension; i++) {
        
        if(p1[i] - p2[i] > epsilon) {
            return false;
        }
        
        if(abs(p1[i] - p2[i]) > epsilon) {
            equal = false;
        }
    }
    
    return !equal;
}


#endif
