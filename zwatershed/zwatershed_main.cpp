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
#include "zwatershed.h"
#include "zwatershed_util/agglomeration.hpp"
#include "zwatershed_util/region_graph.hpp"
#include "zwatershed_util/basic_watershed.hpp"
#include "zwatershed_util/limit_functions.hpp"
#include "zwatershed_util/types.hpp"
#include "zwatershed_util/main_helper.hpp"

// arb funcs
#include "zwatershed_util/region_graph_arb.hpp"
#include "zwatershed_util/basic_watershed_arb.hpp"
#include "zwatershed_util/main_helper_arb.hpp"

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

	// wrap affinities (no copy)
	affinity_graph_ref_ptr<float> affinities(
			new affinity_graph_ref<float>(
					affinity_data,
					boost::extents[width][height][depth][3],
					boost::fortran_storage_order()
			)
	);

	// wrap ground-truth (no copy)
	volume_ref_ptr<uint32_t> ground_truth;
	if (ground_truth_data != 0) {

		ground_truth = volume_ref_ptr<uint32_t>(
				new volume_ref<uint32_t>(
						ground_truth_data,
						boost::extents[width][height][depth],
						boost::fortran_storage_order()
				)
		);
	}

	std::cout << "performing initial watershed segmentation..." << std::endl;

	volume_ptr<uint64_t> segmentation;
	std::vector<size_t> counts;
	std::tie(segmentation, counts) = watershed<uint64_t>(affinities, LOW, HIGH);

	std::cout << "extracting region graph..." << std::endl;
	region_graph_ptr<uint64_t, float> reg_graph = get_region_graph(affinities, segmentation, counts.size() - 1);

	std::vector<Metrics> threshold_metrics;

	for (int i = 0; i < thresholds.size(); i++) {

		size_t threshold = thresholds[i];

		std::cout << "merging until threshold " << threshold << std::endl;
		merge_segments_with_function(
				segmentation,
				reg_graph,
				counts,
				square((double)threshold),
				10, // TODO: what is it?
				RECREATE_RG
		);

		// make a copy of the current segmentation
		std::copy(segmentation->data(), segmentation->data() + num_voxels, &segmentation_data[i][0]);

		if (ground_truth) {

			std::cout << "evaluating current segmentation against ground-truth" << std::endl;

			auto m = compare_volumes(*ground_truth, *segmentation, width, height, depth);
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

ZwatershedResult zwshed_initial_c(const size_t dimX, const size_t dimY, const size_t dimZ, float* affs)
{
    std::cout << "calculating basic watershed..." << std::endl;

    // read data
    volume_ptr<uint64_t> seg_ref;
    std::vector<std::size_t> counts_ref;
    affinity_graph_ptr<float> aff(new affinity_graph<float>
                              (boost::extents[dimX][dimY][dimZ][3],
                               boost::fortran_storage_order()));
    for(size_t i=0;i<dimX*dimY*dimZ*3;i++)
        aff->data()[i] = affs[i];
    std::tie(seg_ref , counts_ref) = watershed<uint64_t>(aff, LOW, HIGH);


    // calculate region graph
    std::cout << "calculating rgn graph..." << std::endl;
    auto rg = get_region_graph(aff, seg_ref , counts_ref.size()-1);

    // save and return
    ZwatershedResult returnMap;
    RegionGraph rg_data;
    for ( const auto& e: *rg ){
        RegionGraphEdge edge = {std::get<1>(e), std::get<2>(e)};
        float weight = std::get<0>(e);
        rg_data.edges.push_back(edge);
        rg_data.weights.push_back(weight);
    }
    std::vector<uint64_t> seg_data;
    std::vector<size_t> counts_data;
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_data.push_back(seg_ref->data()[i]);
    for (const auto& x:counts_ref)
        counts_data.push_back(x);
    returnMap.rg=rg_data;
    returnMap.seg=seg_data;
    returnMap.counts=counts_data;
    return returnMap;
 }


ZwatershedResult merge_with_stats(size_t dimX, size_t dimY, size_t dimZ, uint64_t * gt, RegionGraph& rgn_graph, const std::vector<uint64_t>& seg_in, uint64_t*counts_in, size_t counts_len, size_t thresh){

    size_t rgn_graph_len = rgn_graph.edges.size();

    //read data
    volume_ptr<uint64_t> gt_ptr(new volume<uint64_t> (boost::extents[dimX][dimY][dimZ], boost::fortran_storage_order()));
    volume_ptr<uint64_t> seg(new volume<uint64_t> (boost::extents[dimX][dimY][dimZ], boost::fortran_storage_order()));
    std::vector<std::size_t> counts = * new std::vector<std::size_t>();
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
    for(size_t i=0;i<dimX*dimY*dimZ;i++){
        gt_ptr->data()[i] = gt[i];
        seg->data()[i] = seg_in[i];
    }
    for(size_t i=0;i<counts_len;i++)
        counts.push_back(counts_in[i]);
    for(size_t i=0;i<rgn_graph_len;i++) {

		float weight = rgn_graph.weights[i];
		RegionGraphEdge edge = rgn_graph.edges[i];
        (*rg).emplace_back(weight, edge.id1, edge.id2);
    }

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(seg, rg, counts, square(t), 10,RECREATE_RG);

    // save
	ZwatershedResult returnMap;
    std::vector<uint64_t> seg_vector;
    std::vector<double> stats;
    RegionGraph rg_data;
    std::vector<size_t> counts_data;
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_vector.push_back(seg->data()[i]);
    auto x = compare_volumes(*gt_ptr, *seg, dimX,dimY,dimZ);
    stats.push_back(std::get<0>(x));
    stats.push_back(std::get<1>(x));
    stats.push_back(std::get<2>(x));
    stats.push_back(std::get<3>(x));
    for ( const auto& e: *rg ){
        RegionGraphEdge edge = {std::get<1>(e), std::get<2>(e)};
        float weight = std::get<0>(e);
        rg_data.edges.push_back(edge);
        rg_data.weights.push_back(weight);
    }
    for (const auto& x:counts)
        counts_data.push_back(x);
    returnMap.seg = seg_vector;
    returnMap.stats = stats;
    returnMap.rg=rg_data;
    returnMap.counts = counts_data;
    return returnMap;
}

ZwatershedResult merge_no_stats(size_t dimX, size_t dimY, size_t dimZ, RegionGraph& rgn_graph, const std::vector<uint64_t>& seg_in, uint64_t*counts_in, size_t counts_len, size_t thresh){
    std::cout << "evaluating..." << std::endl;

    size_t rgn_graph_len = rgn_graph.edges.size();

    // read data
    volume_ptr<uint64_t> seg(new volume<uint64_t> (boost::extents[dimX][dimY][dimZ], boost::fortran_storage_order()));
    std::vector<std::size_t> counts = * new std::vector<std::size_t>();
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg->data()[i] = seg_in[i];
    for(size_t i=0;i<counts_len;i++)
        counts.push_back(counts_in[i]);
    for(size_t i=0;i<rgn_graph_len;i++) {

		float weight = rgn_graph.weights[i];
		RegionGraphEdge edge = rgn_graph.edges[i];
        (*rg).emplace_back(weight, edge.id1, edge.id2);
    }

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(seg, rg, counts, square(t), 10,RECREATE_RG);

	// save
	ZwatershedResult returnMap;
    std::vector<uint64_t> seg_vector;
    RegionGraph rg_data;
    std::vector<size_t> counts_data;
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_vector.push_back(seg->data()[i]);
    for ( const auto& e: *rg ){
        RegionGraphEdge edge = {std::get<1>(e), std::get<2>(e)};
        float weight = std::get<0>(e);
        rg_data.edges.push_back(edge);
        rg_data.weights.push_back(weight);
    }
    for (const auto& x:counts)
        counts_data.push_back(x);
    returnMap.seg = seg_vector;
    returnMap.rg=rg_data;
    returnMap.counts = counts_data;
    return returnMap;
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
        rg_data.push_back(std::get<1>(e));
        rg_data.push_back(std::get<2>(e));
        rg_data.push_back(std::get<0>(e));
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
    std::vector<std::size_t> counts = * new std::vector<std::size_t>();
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
    for(size_t i=0;i<dimX*dimY*dimZ;i++){
        gt_ptr->data()[i] = gt[i];
        seg->data()[i] = seg_in[i];
    }
    for(size_t i=0;i<counts_len;i++)
        counts.push_back(counts_in[i]);
    for(size_t i=0;i<rgn_graph_len;i++)
        (*rg).emplace_back(rgn_graph[i*3+2],rgn_graph[i*3],rgn_graph[i*3+1]);

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(seg, rg, counts, square(t), 10,RECREATE_RG);

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
        rg_data.push_back(std::get<1>(e));
        rg_data.push_back(std::get<2>(e));
        rg_data.push_back(std::get<0>(e));
    }
    for (const auto& x:counts)
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
    std::vector<std::size_t> counts = * new std::vector<std::size_t>();
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
    for(size_t i=0;i<dimX*dimY*dimZ;i++){
        seg->data()[i] = seg_in[i];
    }
    for(size_t i=0;i<counts_len;i++)
        counts.push_back(counts_in[i]);
    for(size_t i=0;i<rgn_graph_len;i++)
        (*rg).emplace_back(rgn_graph[i*3+2],rgn_graph[i*3],rgn_graph[i*3+1]);

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(seg, rg, counts, square(t), 10,RECREATE_RG);

    // save
    std::map<std::string,std::vector<double>> returnMap;
    std::vector<double> seg_vector;
    std::vector<double> rg_data; // = * (new std::list<float>());
    std::vector<double> counts_data; // = * (new std::list<float>());
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_vector.push_back(((double)(seg->data()[i])));

    for ( const auto& e: *rg ){
        rg_data.push_back(std::get<1>(e));
        rg_data.push_back(std::get<2>(e));
        rg_data.push_back(std::get<0>(e));
    }
    for (const auto& x:counts)
        counts_data.push_back(x);
    returnMap["seg"] = seg_vector;
    returnMap["rg"]=rg_data;
    returnMap["counts"] = counts_data;
    return returnMap;
}
