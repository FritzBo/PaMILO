//
//  hungarian.cpp
//  mco
//
//  Created by Fritz Bökler on 03.04.14.
//
//

#include <mco/ap/basic/lex_hungarian.h>

#include <functional>
#include <list>
#include <set>
#include <deque>

using std::function;
using std::list;
using std::set;
using std::deque;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::edge;
using ogdf::node;
using ogdf::NodeArray;
using ogdf::EdgeArray;

#include <mco/basic/point.h>

namespace mco {
    
inline void augment(node,
                    NodeArray<node>& exposed_,
                    NodeArray<node>& mate_,
                    NodeArray<node>& label_);
    
inline bool check_equality_subgraph(const Point& cost,
                                    const Point& dual1,
                                    const Point& dual2);
    
inline bool slack_at_least_zero(const Point& cost,
                                const Point& dual1,
                                const Point& dual2);
    
inline bool better_slack(const Point& slack,
                         const Point& cost,
                         const Point& dual1,
                         const Point& dual2);

    

Point LexHungarianMethod::
Solve(const Graph& graph,
      function<Point* (edge)> edge_costs,
      unsigned dimension,
      const set<node>& agents) {
    
    LexPointComparator lex_le(1E-9);
    EqualityPointComparator eq(1E-9);
    
    Point                   null_(dimension);
    Point                   infinity_(numeric_limits<double>::infinity(),
                                      dimension);

    
    NodeArray<list<node>>   A_(graph);
    NodeArray<node>         mate_(graph);
    NodeArray<node>         exposed_(graph);
    NodeArray<node>         neighbour_(graph);
    NodeArray<node>         label_(graph);
    
    NodeArray<Point>        dual_variables_(graph, null_);
    
    NodeArray<Point>        slack_(graph, infinity_);
    
    NodeArray<unsigned>     count_(graph, 0);
    
    deque<node>             queue_;
    
    list<node> non_agents;
    
    edge e;
    node n;
    
    forall_nodes(n, graph) {
        mate_[n] = n;
        if(agents.count(n) == 0) {
            non_agents.push_back(n);
            edge minimum = n->firstAdj()->theEdge();
            forall_adj_edges(e, n) {
                if(lex_le(edge_costs(e), edge_costs(minimum))) {
                    minimum = e;
                }
            }
            
            dual_variables_[n] = *edge_costs(minimum);
        }
    }
    
    bool endstage = false;
    for(unsigned int s = 0; s < agents.size(); ++s) {
        endstage = false;
        
        // Initialize agents
        for(auto agent : agents) {
            A_[agent].clear();
            exposed_[agent] = agent;
            label_[agent] = agent;
            count_[agent] = 0;
        }
        
        // Initialize jobs
        for(auto job : non_agents) {
            slack_[job] = infinity_;
            neighbour_[job] = job;
        }
        
        // Initialize equality subgraph
        forall_edges(e, graph)
        if(check_equality_subgraph(*edge_costs(e),
                                   dual_variables_[e->source()],
                                   dual_variables_[e->target()])) {
            
            auto job = e->target();
            auto agent = e->source();
            if(mate_[job] == job)
            exposed_[agent] = job;
            else if(agent != mate_[job])
            A_[agent].push_back(mate_[job]);
        }
        
        queue_.clear();
        for(auto agent : agents)
        if(mate_[agent] == agent) {
            if(exposed_[agent] != agent) {
                endstage = true;
                augment(agent, exposed_, mate_, label_);
                break;
            }
            queue_.push_back(agent);
            label_[agent] = agent;
            forall_adj_edges(e, agent) {
                if(slack_at_least_zero(*edge_costs(e),
                                       dual_variables_[e->source()],
                                       dual_variables_[e->target()]) &&
                   
                   better_slack(slack_[e->target()],
                                *edge_costs(e),
                                dual_variables_[e->source()],
                                dual_variables_[e->target()])) {
                       
                       slack_[e->target()] = *edge_costs(e) - dual_variables_[e->source()] - dual_variables_[e->target()];
                       if(neighbour_[e->target()] != e->target())
                       count_[neighbour_[e->target()]] -= 1;
                       count_[agent] += 1;
                       neighbour_[e->target()] = agent;
                   }
            }
        }
        
        if(endstage) {
            continue;
        }
        
        while(true) {
            
            while(!queue_.empty()) {
                auto agent = queue_.front();
                queue_.pop_front();
                for(auto next_agent : A_[agent])
                if(label_[next_agent] == next_agent) {
                    label_[next_agent] = agent;
                    
                    if(exposed_[next_agent] != next_agent) {
                        augment(next_agent, exposed_, mate_, label_);
                        endstage = true;
                        break;
                    }
                    
                    queue_.push_back(next_agent);
                    
                    forall_adj_edges(e, next_agent) {
                        if( slack_at_least_zero(*edge_costs(e),
                                                dual_variables_[e->source()],
                                                dual_variables_[e->target()]) &&
                           better_slack(slack_[e->target()],
                                        *edge_costs(e),
                                        dual_variables_[e->source()],
                                        dual_variables_[e->target()])) {
                               
                               slack_[e->target()] = *edge_costs(e) - dual_variables_[e->source()] - dual_variables_[e->target()];
                               if(neighbour_[e->target()] != e->target())
                               count_[neighbour_[e->target()]] -= 1;
                               count_[next_agent] += 1;
                               neighbour_[e->target()] = next_agent;
                           }
                    }
                }
                
                if(endstage) {
                    break;
                }
            }
            
            if(endstage)
            break;
            
            Point minimum_slack = infinity_;
            for(auto job : non_agents) {
                if(lex_le(null_, slack_[job]) &&
                   lex_le(slack_[job], minimum_slack)) {
                    
                    minimum_slack = slack_[job];
                }
            }
            Point theta = minimum_slack;
            theta *= 0.5;
            
            for(auto agent : agents) {
                if(label_[agent] != agent || count_[agent] > 0) {
                    dual_variables_[agent] += theta;
                } else {
                    dual_variables_[agent] -= theta;
                }
            }
            
            for(auto job : non_agents) {
                if(eq(slack_[job], null_)) {
                    dual_variables_[job] -= theta;
                } else {
                    dual_variables_[job] += theta;
                }
            }
            
            theta *= 2;
            for(auto job : non_agents)
            if(lex_le(null_, slack_[job])) {
                slack_[job] -= theta;
                if(eq(slack_[job], null_)) {
                    if(mate_[job] == job) {
                        exposed_[neighbour_[job]] = job;
                        augment(neighbour_[job], exposed_, mate_, label_);
                        endstage = true;
                        break;
                    } else {
                        queue_.push_back(neighbour_[job]);
                        A_[neighbour_[job]].push_back(mate_[job]);
                    }
                }
            }
            
            if(endstage) {
                break;
            }
        }
    }
    
#ifndef NDEBUG
//    for(auto n : graph.nodes) {
//        if(n->index() < agents.size()) {
//            std::cout << "(" << n << ", " << mate_[n] << "), ";
//        }
//    }
//    std::cout << "-- ";
#endif
    
    if(value_.dimension() != dimension) {
        value_ = Point(dimension);
    } else {
        for(unsigned i = 0; i < dimension; ++i) {
            value_[i] = 0.0;
        }
    }
    
    for(auto n : graph.nodes) {
        value_ += dual_variables_[n];
    }
    
#ifndef NDEBUG
//    std::cout << value_ << std::endl;
#endif
    
    return value_;
}

inline bool check_equality_subgraph(const Point &cost,
                                    const Point &dual1,
                                    const Point &dual2) {
    
    for(unsigned i = 0; i < cost.dimension(); ++i) {
        if(abs(cost[i] - (dual1[i] + dual2[i])) > 1E-9) {
            return false;
        }
    }
    
    return true;
}

inline bool slack_at_least_zero(const Point &cost,
                                const Point &dual1,
                                const Point &dual2) {
    
    for(unsigned int i = 0; i < cost.dimension(); ++i) {
        if(cost[i] - dual1[i] - dual2[i] > 1E-9) {
            return true;
        } else if(cost[i] - dual1[i] - dual2[i] < -1E-9) {
            return false;
        }
    }
    
    return true;
}

inline bool better_slack(const Point &slack,
                         const Point &cost,
                         const Point &dual1,
                         const Point &dual2) {
    
    for(unsigned int i = 0; i < cost.dimension(); ++i) {
        if(slack[i] - (cost[i] - dual1[i] - dual2[i]) > 1E-9) {
            return true;
        } else if(slack[i] - (cost[i] - dual1[i] - dual2[i]) < -1E-9) {
            return false;
        }
    }
    return false;
}

inline void augment(node v,
                    NodeArray<node>& exposed_,
                    NodeArray<node>& mate_,
                    NodeArray<node>& label_) {
    
    while(label_[v] != v) {
        exposed_[label_[v]] = mate_[v];
        mate_[v] = exposed_[v];
        mate_[exposed_[v]] = v;
        v = label_[v];
    }
    mate_[v] = exposed_[v];
    mate_[exposed_[v]] = v;
}

}