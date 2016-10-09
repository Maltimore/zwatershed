#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <vector>
#include <utility>

struct Metrics {

	double voi_split;
	double voi_merge;
	double rand_split;
	double rand_merge;
};

std::vector<Metrics> process_thresholds(
		const std::vector<size_t>& thresholds,
		size_t width, size_t height, size_t depth,
		const float* affinity_data,
		const std::vector<uint64_t*>& segmentation_data,
		const uint32_t* ground_truth_data);


struct RegionGraphEdge {

	uint64_t id1;
	uint64_t id2;
};

struct RegionGraph {

	std::vector<RegionGraphEdge> edges;
	std::vector<float>           weights;
};

struct ZwatershedResult {

	RegionGraph rg;
	std::vector<uint64_t> seg;
	std::vector<size_t> counts;
	std::vector<double> stats;
};

ZwatershedResult zwshed_initial_c(const size_t dx, const size_t dy, const size_t dz, float* affs);

ZwatershedResult merge_with_stats(size_t dx,size_t dy, size_t dz, uint64_t * gt, RegionGraph& rgn_graph, const std::vector<uint64_t>& seg_in, uint64_t*counts, size_t counts_len, size_t thresh);

ZwatershedResult merge_no_stats(size_t dimX, size_t dimY, size_t dimZ, RegionGraph& rgn_graph, const std::vector<uint64_t>& seg_in, uint64_t*counts, size_t counts_len, size_t thresh);

std::map<std::string,std::list<float>> zwshed_initial_c_arb(const size_t dx, const size_t dy, const size_t dz, const uint64_t*node1,
                                               const uint64_t*node2, const float*edgeWeight, const size_t n_edge);

std::map<std::string,std::vector<double>> merge_with_stats_arb(size_t dx,size_t dy, size_t dz, uint64_t * gt, float * rgn_graph,
                                        size_t rgn_graph_len, uint64_t * seg_in, uint64_t*counts, size_t counts_len, size_t thresh);

std::map<std::string,std::vector<double>> merge_no_stats_arb(size_t dx,size_t dy, size_t dz, float * rgn_graph,
                                        size_t rgn_graph_len, uint64_t * seg_in, uint64_t*counts, size_t counts_len, size_t thresh);

#endif
