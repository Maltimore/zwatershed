#ifndef ZWATERSHED_H
#define ZWATERSHED_H

#include <list>
#include <string>
#include <map>
#include <vector>

#include "backend/types.hpp"

struct Metrics {

	double voi_split;
	double voi_merge;
	double rand_split;
	double rand_merge;
};

struct ZwatershedState {

	volume_ref_ptr<uint64_t> segmentation;
	counts_ptr<size_t> counts;
	std::shared_ptr<RegionGraph<uint64_t>> region_graph;
	std::shared_ptr<typename RegionGraph<uint64_t>::template EdgeMap<float>> edge_affinities;
};

std::vector<Metrics> process_thresholds(
		const std::vector<float>& thresholds,
		size_t width, size_t height, size_t depth,
		const float* affinity_data,
		const std::vector<uint64_t*>& segmentation_data,
		const uint32_t* ground_truth_data = 0);

ZwatershedState get_initial_state(
		size_t width, size_t height, size_t depth,
		const float* affinity_data,
		uint64_t* segmentation_data);

#endif
