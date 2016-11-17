/* Connected components
 * developed and maintained by Srinivas C. Turaga <sturaga@mit.edu>
 * do not distribute without permission.
 */

#include <memory>
#include <type_traits>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <chrono>
#include <string>

//#pragma once
#include "c_frontend.h"
#include "backend/agglomeration.hpp"
#include "backend/region_graph.hpp"
#include "backend/basic_watershed.hpp"
#include "backend/limit_functions.hpp"
#include "backend/main_helper.hpp"
#include "backend/IterativeRegionMerging.hpp"
#include "backend/MergeFunctions.hpp"

using namespace std;
// these values based on 5% at iter = 10000
double LOW=  .0001;
double HIGH= .9999;
bool RECREATE_RG = true;

std::vector<Metrics> process_thresholds(
		const std::vector<float>& thresholds,
		size_t width, size_t height, size_t depth,
		const float* affinity_data,
		const std::vector<uint64_t*>& segmentation_data,
		const uint32_t* ground_truth_data) {

	size_t num_voxels = width*height*depth;

	assert(thresholds.size() == segmentation_data.size());

	ZwatershedState state = get_initial_state(
			width, height, depth,
			affinity_data,
			segmentation_data[0]);

	std::vector<Metrics> threshold_metrics;

	IterativeRegionMerging<uint64_t, float> regionMerging(state.region_graph);
	MedianAffinity<uint64_t, float> mergeFunction(*state.edge_affinities);

	for (int i = 0; i < thresholds.size(); i++) {

		float threshold = thresholds[i];

		std::cout << "merging until threshold " << threshold << std::endl;
		regionMerging.mergeUntil(
				mergeFunction,
				threshold);

		std::cout << "extracting segmentation" << std::endl;

		// wrap segmentation for current iteration (no copy)
		volume_ref<uint64_t> current_segmentation(
				segmentation_data[i],
				boost::extents[width][height][depth]
		);

		regionMerging.extractSegmentation(current_segmentation);

		// make a copy of the current segmentation for the next iteration
		if (i < segmentation_data.size() - 1)
			std::copy(segmentation_data[i], segmentation_data[i] + num_voxels, segmentation_data[i+1]);

		if (ground_truth_data != 0) {

			std::cout << "evaluating current segmentation against ground-truth" << std::endl;

			// wrap ground-truth (no copy)
			volume_const_ref<uint32_t> ground_truth(
					ground_truth_data,
					boost::extents[width][height][depth]
			);

			auto m = compare_volumes(ground_truth, current_segmentation, width, height, depth);
			Metrics metrics;
			metrics.rand_split = std::get<0>(m);
			metrics.rand_merge = std::get<1>(m);
			metrics.voi_split  = std::get<2>(m);
			metrics.voi_merge  = std::get<3>(m);

			threshold_metrics.push_back(metrics);
		}
	}

	return threshold_metrics;
}

ZwatershedState get_initial_state(
		size_t width, size_t height, size_t depth,
		const float* affinity_data,
		uint64_t* segmentation_data) {

	size_t num_voxels = width*height*depth;

	// wrap affinities (no copy)
	affinity_graph_ref<float> affinities(
			affinity_data,
			boost::extents[3][width][height][depth]
	);

	// wrap segmentation array (no copy)
	volume_ref_ptr<uint64_t> segmentation(
			new volume_ref<uint64_t>(
					segmentation_data,
					boost::extents[width][height][depth]
			)
	);

	// create counts data structure
	counts_ptr<size_t> counts(new counts_t<size_t>());

	std::cout << "performing initial watershed segmentation..." << std::endl;

	watershed(affinities, LOW, HIGH, *segmentation, *counts);

	std::cout << "extracting region graph..." << std::endl;

	std::size_t numNodes = counts->size();

	std::shared_ptr<RegionGraph<uint64_t>> regionGraph(
			new RegionGraph<uint64_t>(numNodes)
	);

	std::shared_ptr<RegionGraph<uint64_t>::EdgeMap<float>> edgeAffinities(
			new RegionGraph<uint64_t>::EdgeMap<float>(*regionGraph)
	);

	get_region_graph(
			affinities,
			*segmentation,
			counts->size() - 1,
			*regionGraph,
			*edgeAffinities);

	ZwatershedState initial_state;
	initial_state.region_graph = regionGraph;
	initial_state.edge_affinities = edgeAffinities;
	initial_state.segmentation = segmentation;
	initial_state.counts = counts;

	return initial_state;
}

