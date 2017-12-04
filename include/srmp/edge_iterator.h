#ifndef SRMP_EDGE_ITERATOR_H
#define SRMP_EDGE_ITERATOR_H

#include <functional>
#include <vector>

#include <srmp/SRMP.h>

namespace srmpLib {

class EdgeIterator {
public:
	struct Edge {
		std::vector<int> alpha;
		std::vector<int> beta;
		double *message;
		size_t message_size;

		double* message_begin() const { return message; }
		double* message_end() const { return message + message_size; }
	};

	EdgeIterator(const Energy *energy)
	: energy_(energy)
	{
		edge_ = energy_->edges->ScanFirst(iterator_);
	}

	void operator++()
	{
		edge_ = energy_->edges->ScanNext(iterator_);
	}

	Edge* operator->()
	{
		convert_factor(edge_->A, &current_.alpha);
		convert_factor(edge_->B, &current_.beta);
		current_.message = edge_->m;
		current_.message_size = edge_->B->K;
		return &current_;
	}

	Edge& operator*()
	{
		return *this->operator->();
	}

	bool valid() const
	{
		return edge_ != NULL;
	}

private:
	const Energy *energy_;
	const Energy::Edge * edge_;
	Block<Energy::Edge>::iterator iterator_;
	Edge current_;

	void convert_factor(const Energy::Factor *factor, std::vector<int> *out) const
	{
		// Note that the differentation of Factor and NonSingletonFactor hits
		// us here, but this functions seems to handle both cases correctly.
		// Unfortunately it does not obey const correctness... Also the
		// function takes a REFERENCE TO A POINTER, so we have to store it in
		// an rvalue. :(
		Energy::Factor *tmp = const_cast<Energy::Factor*>(factor);

		Energy::Node **nodes = energy_->GetSortedNodesPtr(tmp);
		out->resize(factor->arity);
		for (int i = 0; i < factor->arity; ++i)
			(*out)[i] = nodes[i] - energy_->nodes;
	}
};

}

#endif
