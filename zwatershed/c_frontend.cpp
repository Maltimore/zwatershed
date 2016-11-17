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

// arb funcs
#include "backend/region_graph_arb.hpp"
#include "backend/basic_watershed_arb.hpp"
#include "backend/main_helper_arb.hpp"

using namespace std;
// these values based on 5% at iter = 10000
double LOW=  .0001;
double HIGH= .9999;
bool RECREATE_RG = true;

std::vector<Metrics> process_thresholds(
		const std::vector<size_t>& thresholds,
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
	MedianAffinity mergeFunction;

	for (int i = 0; i < thresholds.size(); i++) {

		size_t threshold = thresholds[i];

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
	region_graph_ptr<uint64_t, float> reg_graph(new region_graph<uint64_t, float>());
	get_region_graph(affinities, *segmentation, counts->size() - 1, *reg_graph);

	ZwatershedState initial_state;
	initial_state.region_graph = reg_graph;
	initial_state.segmentation = segmentation;
	initial_state.counts = counts;

	return initial_state;
}


/////////////////////////////////////////arb nhoods/////////////////////////////////////////



std::map<std::string,std::list<float>> zwshed_initial_c_arb(const size_t dimX, const size_t dimY, const size_t dimZ, const uint64_t*node1,
                                               const uint64_t*node2, const float*edgeWeight, const size_t n_edge){
    // read data
    std::cout << "calculating basic watershed..." << std::endl;
    volume_ptr<uint64_t> seg_ref;
    std::vector<std::size_t> counts_ref;
    std::tie(seg_ref , counts_ref) = watershed_arb<uint64_t>(dimX,dimY,dimZ,node1, node2, edgeWeight, n_edge, LOW, HIGH);
    auto seg = *seg_ref;

    // calculate region graph
    std::cout << "calculating rgn graph..." << std::endl;
    auto rg = get_region_graph_arb(node1, node2, edgeWeight, n_edge, seg_ref , counts_ref.size()-1);

    // save and return
    std::map<std::string,std::list<float>> returnMap;

    std::list<float> rg_data = * (new std::list<float>());
    for ( const auto& e: *rg ){
		// TODO: unsafe cast from size_t to float
        rg_data.push_back(e.id1);
        rg_data.push_back(e.id2);
        rg_data.push_back(e.affinity);
    }
    std::list<float> seg_data = * (new std::list<float>());
    std::list<float> counts_data = * (new std::list<float>());
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_data.push_back(seg_ref->data()[i]);
    for (const auto& x:counts_ref)
        counts_data.push_back(x);

    returnMap["rg"]=rg_data;
    returnMap["seg"]=seg_data;
    returnMap["counts"]=counts_data;

    return returnMap;
 }


std::map<std::string,std::vector<double>> merge_with_stats_arb(size_t dimX,size_t dimY, size_t dimZ, uint64_t * gt, float * rgn_graph,
size_t rgn_graph_len, uint64_t * seg_in, uint64_t*counts_in, size_t counts_len, size_t thresh){

    //read data
    volume_ptr<uint64_t> gt_ptr(new volume<uint64_t> (boost::extents[dimX][dimY][dimZ], boost::c_storage_order() )); //, boost::fortran_storage_order()));
    volume_ptr<uint64_t> seg(new volume<uint64_t> (boost::extents[dimX][dimY][dimZ], boost::c_storage_order()));
    counts_ptr<std::size_t> counts(new counts_t<std::size_t>());
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
    for(size_t i=0;i<dimX*dimY*dimZ;i++){
        gt_ptr->data()[i] = gt[i];
        seg->data()[i] = seg_in[i];
    }
    for(size_t i=0;i<counts_len;i++)
        counts->push_back(counts_in[i]);
    for(size_t i=0;i<rgn_graph_len;i++)
        (*rg).emplace_back(rgn_graph[i*3+2],rgn_graph[i*3],rgn_graph[i*3+1]);

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(*seg, rg, *counts, square(t), 10,RECREATE_RG);

    // save
    std::map<std::string,std::vector<double>> returnMap;
    std::vector<double> seg_vector;
    std::vector<double> r;
    std::vector<double> rg_data; // = * (new std::list<float>());
    std::vector<double> counts_data; // = * (new std::list<float>());
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_vector.push_back(((double)(seg->data()[i])));
    auto x = compare_volumes_arb(*gt_ptr, *seg, dimX,dimY,dimZ);
    r.push_back(x.first);
    r.push_back(x.second);
    for ( const auto& e: *rg ){
		// TODO: unsafe cast from size_t to float
        rg_data.push_back(e.id1);
        rg_data.push_back(e.id2);
        rg_data.push_back(e.affinity);
    }
    for (const auto& x:*counts)
        counts_data.push_back(x);
    returnMap["seg"] = seg_vector;
    returnMap["stats"] = r;
    returnMap["rg"]=rg_data;
    returnMap["counts"] = counts_data;

    return returnMap;
}

std::map<std::string,std::vector<double>> merge_no_stats_arb(size_t dimX,size_t dimY, size_t dimZ, float * rgn_graph,
size_t rgn_graph_len, uint64_t * seg_in, uint64_t*counts_in, size_t counts_len, size_t thresh){

    //read data
    volume_ptr<uint64_t> seg(new volume<uint64_t> (boost::extents[dimX][dimY][dimZ]));
    counts_ptr<std::size_t> counts(new counts_t<std::size_t>());
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
    for(size_t i=0;i<dimX*dimY*dimZ;i++){
        seg->data()[i] = seg_in[i];
    }
    for(size_t i=0;i<counts_len;i++)
        counts->push_back(counts_in[i]);
    for(size_t i=0;i<rgn_graph_len;i++)
        (*rg).emplace_back(rgn_graph[i*3+2],rgn_graph[i*3],rgn_graph[i*3+1]);

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(*seg, rg, *counts, square(t), 10,RECREATE_RG);

    // save
    std::map<std::string,std::vector<double>> returnMap;
    std::vector<double> seg_vector;
    std::vector<double> rg_data; // = * (new std::list<float>());
    std::vector<double> counts_data; // = * (new std::list<float>());
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_vector.push_back(((double)(seg->data()[i])));

    for ( const auto& e: *rg ){
		// TODO: unsafe cast from size_t to float
        rg_data.push_back(e.id1);
        rg_data.push_back(e.id2);
        rg_data.push_back(e.affinity);
    }
    for (const auto& x:*counts)
        counts_data.push_back(x);
    returnMap["seg"] = seg_vector;
    returnMap["rg"]=rg_data;
    returnMap["counts"] = counts_data;
    return returnMap;
}
