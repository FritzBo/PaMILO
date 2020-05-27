//
//  utility.h
//  mco
//
//  Created by Fritz Bökler on 09.04.14.
//
//

#ifndef __mco__utility__
#define __mco__utility__

namespace mco {

class NonDominatedPartition {
public:
    NonDominatedPartition(double epsilon = 1E-8)
    :   epsilon_(epsilon) {}
    
    template<typename ConstIterator, typename BackInserter>
    bool operator()(ConstIterator source1_begin,
                    ConstIterator source1_end,
                    ConstIterator source2_begin,
                    ConstIterator source2_end,
                    BackInserter nondominated_subset,
                    BackInserter dominated_subset);
private:
    double epsilon_;
};

template<typename ConstIterator, typename BackInserter>
bool NonDominatedPartition::
operator()(ConstIterator source1_begin,
                ConstIterator source1_end,
                ConstIterator source2_begin,
                ConstIterator source2_end,
                BackInserter nondominated_subset,
                BackInserter dominated_subset) {
    
	bool new_labels = false;
    
    //	if(dim == 2) {
    //
    //		auto it_s1 = source1.cbegin();
    //		auto it_s2 = source2.cbegin();
    //
    //		if((**it_s1)[0] <= (**it_s2)[0] && (**it_s1)[1] > (**it_s2)[1]) {
    //			nondominated_subset.push_back(*it_s1);
    //			new_labels = true;
    //			it_s1++;
    //		} else {
    //			nondominated_subset.push_back(*it_s2);
    //			it_s2++;
    //		}
    //
    //		while(!(it_s1 == source1.end() && it_s2 == source2.end())) {
    //			if(it_s2 != source2.end() && (it_s1 == source1.end() || ((**it_s1)[0] >= (**it_s2)[0] && ((**it_s1)[1] < (**it_s2)[1])))) {
    //				if((**it_s2)[1] <= (*nondominated_subset.back())[1] && (**it_s2) != (*nondominated_subset.back()))
    //					nondominated_subset.push_back(*it_s2);
    //				else
    //					dominated_subset.push_back(*it_s2);
    //				it_s2++;
    //			} else {
    //				if((**it_s1)[1] <= (*nondominated_subset.back())[1] && (**it_s1) != (*nondominated_subset.back())) {
    //					nondominated_subset.push_back(*it_s1);
    //					new_labels = true;
    //				} else
    //					dominated_subset.push_back(*it_s1);
    //				it_s1++;
    //			}
    //		}
    //
    //	} else {
    
    using Iterator1 = decltype(source1_begin);
    using Iterator2 = decltype(source2_begin);
    
    unsigned source2_size = source2_end - source2_begin;
    
    vector<bool> marker(source2_size);
    
    EqualityPointComparator eq_comp(epsilon_);
    ComponentwisePointComparator comp_leq(epsilon_, false);
    ParetoDominationPointComparator pareto_dominates(epsilon_);
    
    for(unsigned i = 0; i < source2_size; ++i) {
        marker[i] = false;
    }
    
    bool dominated;
    
    for(Iterator1 it1 = source1_begin; it1 < source2_end; ++it1) {
        dominated = false;
        auto label_source1 = *it1;
        
        unsigned int i = 0;
        for(Iterator2 it2 = source2_begin; it2 < source2_end; ++it2) {
            auto label_source2 = *it2;
            
            if(eq_comp(label_source1, label_source2)) {
                dominated = true;
                break;
            }
            
            if(comp_leq(label_source1, label_source2)) {
                marker[i] = true;
            }
            
            if(comp_leq(label_source2, label_source2)) {
                dominated = true;
                break;
            }
            
            i++;
        }
        
        if(!dominated) {
            nondominated_subset.push_back(label_source1);
            new_labels = true;
        } else
            dominated_subset.push_back(label_source1);
    }
    
    unsigned int i = 0;
    
    for(Iterator2 it2 = source2_begin; it2 < source2_end; ++it2) {
        auto label = *it2;
        
        if(!marker[i])
            nondominated_subset.push_back(label);
        else
            dominated_subset.push_back(label);
        i++;
    }
    
	return new_labels;
}

}
    
#endif /* defined(__mco__utility__) */
