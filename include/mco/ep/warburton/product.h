#pragma once
/*
 * product_iterator.h
 *
 *  Created on: 08.04.2013
 *      Author: fritz
 */

#ifndef PRODUCT_H_
#define PRODUCT_H_

#include <utility>
#include <iterator>
#include <vector>

#include <mco/point.h>

namespace mco {

template<typename InputIterator, typename OutputIterator>
class Product {

	std::vector<std::pair<InputIterator, InputIterator>> ranges_;
	OutputIterator begin_;
	std::vector<InputIterator> current_;
	unsigned int number_ranges_;
	bool end_;

public:
	Product(std::vector<std::pair<InputIterator, InputIterator>> ranges, OutputIterator begin) : ranges_(ranges), begin_(begin), end_(false) {
		number_ranges_ = ranges.size();
		current_.reserve(number_ranges_);
		for(unsigned int i = 0; i < number_ranges_; ++i)
			current_[i] = ranges_[i].first;
	}

	void operator()() {

		if(end_ == true)
			return;

		OutputIterator output(begin_);

		for(auto iterator : current_) {
			*output = *iterator;
			output++;
		}

		for(unsigned int i = 0; i < number_ranges_; ++i) {
			if(current_[i] != ranges_[i].second) {
				current_[i]++;
				break;
			} else {
				current_[i] = ranges_[i].first;
				if(i == number_ranges_ - 1)
					end_ = true;
			}
		}
	}

	bool has_next() const {
		return !end_;
	}

};

template<typename InputIterator, typename OutputIterator>
Product<InputIterator, OutputIterator> generate_product(std::vector<std::pair<InputIterator, InputIterator>> ranges, OutputIterator begin) {
	return Product<InputIterator, OutputIterator>(ranges, begin);
}

}

#endif /* PRODUCT_H_ */
