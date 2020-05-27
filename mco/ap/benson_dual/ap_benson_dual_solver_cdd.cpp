/*
 * ap_benson_dual_solver.cpp
 *
 *  Created on: 14.10.2013
 *      Author: fritz
 */

#include <mco/ap/benson_dual/ap_benson_dual_solver_cdd.h>

#include <limits>
#include <cmath>

using std::numeric_limits;
using std::set;
using std::list;
using std::abs;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::edge;
using ogdf::node;
using ogdf::EdgeArray;
using ogdf::NodeArray;

namespace mco {

APBensonDualSolverCDD::APBensonDualSolverCDD(std::shared_ptr<AssignmentInstance> inpt_instance, double epsilon):
	AbstractAPSolver(inpt_instance), epsilon_(epsilon),
	A_(*instance().graph()),
	mate_(*instance().graph()),
	exposed_(*instance().graph()),
	neighbour_(*instance().graph()),
	label_(*instance().graph()),
	dual_variables_(*instance().graph(), Point(0.0, instance().dimension() + 1)),
	slack_(*instance().graph(), Point(numeric_limits<double>::infinity(), instance().dimension() + 1)),
	count_(*instance().graph(), 0),
	edge_costs(*instance().graph(), Point(instance().dimension() + 1)),
	null_(0.0, instance().dimension() + 1),
	infinity_(numeric_limits<double>::infinity(), instance().dimension() + 1),
	cycles_(0) {

	edge e;
	forall_edges(e, *instance().graph())
		for(unsigned int i = 0; i < instance().dimension(); ++i)
			edge_costs[e][i + 1] = (*(*instance().weights())[e])[i];


	dual_benson_solver_ = new DualBensonScalarizerCDD(*this, instance().dimension(), epsilon);
}

APBensonDualSolverCDD::~APBensonDualSolverCDD() {
}

void APBensonDualSolverCDD::augment(node v) {
	while(label_[v] != v) {
		exposed_[label_[v]] = mate_[v];
		mate_[v] = exposed_[v];
		mate_[exposed_[v]] = v;
		v = label_[v];
	}
	mate_[v] = exposed_[v];
	mate_[exposed_[v]] = v;
}

double APBensonDualSolverCDD::Solve_scalarization(Point& weights, Point& value) {

	clock_t start = clock();

	unsigned int dim = instance().dimension();
	Graph &graph = *instance().graph();
	set<node> agents = *instance().agents();
	list<node> non_agents;

	edge e;
	node n;

	forall_edges(e, graph)
		edge_costs[e][0] = *(*instance().weights())[e] * weights;

	forall_nodes(n, graph) {
		mate_[n] = n;
		if(agents.count(n) == 0) {
			non_agents.push_back(n);
			edge minimum = n->firstAdj()->theEdge();
			forall_adj_edges(e, n)
				if(edge_costs[e].is_lexicographic_less(edge_costs[minimum], epsilon_))
					minimum = e;
			dual_variables_[n] = edge_costs[minimum];
		} else {
			dual_variables_[n] = null_;
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
			if(check_equality_subgraph(edge_costs[e], dual_variables_[e->source()], dual_variables_[e->target()], epsilon_)) {
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
					augment(agent);
					break;
				}
				queue_.push_back(agent);
				label_[agent] = agent;
				forall_adj_edges(e, agent) {
					if( slack_at_least_zero(edge_costs[e], dual_variables_[e->source()], dual_variables_[e->target()], epsilon_) &&
						better_slack(slack_[e->target()], edge_costs[e], dual_variables_[e->source()], dual_variables_[e->target()], epsilon_)) {
						slack_[e->target()] = edge_costs[e] - dual_variables_[e->source()] - dual_variables_[e->target()];
						if(neighbour_[e->target()] != e->target())
							count_[neighbour_[e->target()]] -= 1;
						count_[agent] += 1;
						neighbour_[e->target()] = agent;
					}
				}
			}

		if(endstage)
			continue;

		while(true) {

			while(!queue_.empty()) {
				auto agent = queue_.front();
				queue_.pop_front();
				for(auto next_agent : A_[agent])
					if(label_[next_agent] == next_agent) {
						label_[next_agent] = agent;

						if(exposed_[next_agent] != next_agent) {
							augment(next_agent);
							endstage = true;
							break;
						}

						queue_.push_back(next_agent);

						forall_adj_edges(e, next_agent) {
							if( slack_at_least_zero(edge_costs[e], dual_variables_[e->source()], dual_variables_[e->target()], epsilon_) &&
								better_slack(slack_[e->target()], edge_costs[e], dual_variables_[e->source()], dual_variables_[e->target()], epsilon_)) {
								slack_[e->target()] = edge_costs[e] - dual_variables_[e->source()] - dual_variables_[e->target()];
								if(neighbour_[e->target()] != e->target())
									count_[neighbour_[e->target()]] -= 1;
								count_[next_agent] += 1;
								neighbour_[e->target()] = next_agent;
							}
						}
					}

				if(endstage)
					break;
			}

			if(endstage)
				break;

			Point minimum_slack = infinity_;
			for(auto job : non_agents) {
				if(		null_.is_lexicographic_less(slack_[job], epsilon_) &&
						slack_[job].is_lexicographic_less(minimum_slack, epsilon_))
					minimum_slack = slack_[job];
			}
			Point theta = minimum_slack;
			theta *= 0.5;

			for(auto agent : agents) {
				if(label_[agent] != agent || count_[agent] > 0)
					dual_variables_[agent] += theta;
				else
					dual_variables_[agent] -= theta;
			}

			for(auto job : non_agents) {
				if(slack_[job].is_equal(null_, epsilon_))
					dual_variables_[job] -= theta;
				else
					dual_variables_[job] += theta;
			}

			theta *= 2;
			for(auto job : non_agents)
				if(null_.is_lexicographic_less(slack_[job], epsilon_)) {
					slack_[job] -= theta;
					if(slack_[job].is_equal(null_, epsilon_)) {
						if(mate_[job] == job) {
							exposed_[neighbour_[job]] = job;
							augment(neighbour_[job]);
							endstage = true;
							break;
						} else {
							queue_.push_back(neighbour_[job]);
							A_[neighbour_[job]].push_back(mate_[job]);
						}
					}
				}

			if(endstage)
				break;
		}
	}


	forall_nodes(n, graph)
		for(unsigned int i = 0; i < dim; ++i)
			value[i] += dual_variables_[n][i + 1];

	cycles_ += clock() - start;
//	cout << "Solving took " << (clock() - start) / (double) CLOCKS_PER_SEC << "s." << endl;

	return weights * value;
}

bool APBensonDualSolverCDD::check_equality_subgraph(Point &cost, Point &dual1, Point &dual2, double epsilon_) {
	for(unsigned int i = 0; i < cost.dimension(); ++i)
		if(abs(cost[i] - (dual1[i] + dual2[i])) > epsilon_)
			return false;

	return true;
}

bool APBensonDualSolverCDD::slack_at_least_zero(Point &cost, Point &dual1, Point &dual2, double epsilon_) {
	for(unsigned int i = 0; i < cost.dimension(); ++i)
		if(cost[i] - dual1[i] - dual2[i] > epsilon_)
			return true;
		else if(cost[i] - dual1[i] - dual2[i] < -epsilon_)
			return false;

	return true;
}

bool APBensonDualSolverCDD::better_slack(Point &slack, Point &cost, Point &dual1, Point &dual2, double epsilon_) {
	for(unsigned int i = 0; i < cost.dimension(); ++i)
		if(slack[i] - (cost[i] - dual1[i] - dual2[i]) > epsilon_)
			return true;
		else if(slack[i] - (cost[i] - dual1[i] - dual2[i]) < -epsilon_)
			return false;

	return false;
}


} /* namespace mco */
