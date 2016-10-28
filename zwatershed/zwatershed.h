#ifndef ZWATERSHED_H
#define ZWATERSHED_H

#include <list>
#include <string>
#include <map>
#include <vector>

#include "zwatershed_util/types.hpp"

using RegionGraphEdge = region_graph_edge_t<float, uint64_t>;

struct Metrics {

	double voi_split;
	double voi_merge;
	double rand_split;
	double rand_merge;
};

struct ZwatershedState {

	volume_ptr<uint64_t> segmentation;
	counts_ptr<size_t> counts;
	region_graph_ptr<uint64_t, float> region_graph;
};

std::vector<Metrics> process_thresholds(
		const std::vector<size_t>& thresholds,
		size_t width, size_t height, size_t depth,
		const float* affinity_data,
		const std::vector<uint64_t*>& segmentation_data,
		const uint32_t* ground_truth_data = 0);

ZwatershedState get_initial_state(
		size_t width, size_t height, size_t depth,
		const float* affinity_data);

std::map<std::string,std::list<float>> zwshed_initial_c_arb(const size_t dx, const size_t dy, const size_t dz, const uint64_t*node1,
                                               const uint64_t*node2, const float*edgeWeight, const size_t n_edge);

std::map<std::string,std::vector<double>> merge_with_stats_arb(size_t dx,size_t dy, size_t dz, uint64_t * gt, float * rgn_graph,
                                        size_t rgn_graph_len, uint64_t * seg_in, uint64_t*counts, size_t counts_len, size_t thresh);

std::map<std::string,std::vector<double>> merge_no_stats_arb(size_t dx,size_t dy, size_t dz, float * rgn_graph,
                                        size_t rgn_graph_len, uint64_t * seg_in, uint64_t*counts, size_t counts_len, size_t thresh);

#endif
