#pragma once

#include "types.hpp"

#include <cstddef>
#include <iostream>
#include <map>
using namespace std;
template< typename ID, typename F >
inline region_graph_ptr<ID>
get_region_graph_arb( //const affinity_graph_ptr<F>& aff_ptr,
                const ID*node1, const ID*node2, const F*edgeWeight, size_t n_edge,
                  const volume_ptr<ID> seg_ptr, std::size_t max_segid)
{

    region_graph_ptr<ID> rg_ptr( new RegionGraph<ID>() );
	// FIXME: needs to be adapted to new RegionGraph
    //region_graph<ID,F>& rg = *rg_ptr;
    //ID* seg = seg_ptr->data();
    //std::vector<std::map<ID,F>> edges(max_segid+1);
    //std::ptrdiff_t x = 1;
    //for(size_t i=0;i<n_edge;i++){
        //ID n1 = node1[i];
        //ID n2 = node2[i];
        //F w = edgeWeight[i];
        //auto mm = std::minmax(seg[n1],seg[n2]); //need to make sure this is in right order
        //F& curr = edges[mm.first][mm.second];
        //curr = std::max(curr, w);
    //}


    //for ( ID id1 = 1; id1 <= max_segid; ++id1 )
        //for ( const auto& p: edges[id1] )
            //rg.emplace_back(id1, p.first, p.second);


    //std::stable_sort(std::begin(rg), std::end(rg), std::greater<RegionGraphEdge<ID,F>>());


    return rg_ptr;
}

