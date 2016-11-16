#pragma once

#include "types.hpp"

#include <zi/disjoint_sets/disjoint_sets.hpp>
#include <map>
#include <vector>
#include <set>

template<typename V, typename F, typename FN>
inline void
merge_segments_with_function(
        V& seg,
        const region_graph_ptr<typename V::element, F> rg_ptr,
        counts_t<std::size_t>& counts,
        const FN& func,
        size_t low,
        bool recreate_rg)
{
    typedef typename V::element ID;

	// one "set" for each region
	// seems to be the connected components
    zi::disjoint_sets<ID> sets(counts.size());

    region_graph<ID,F>& rg  = *rg_ptr;

	// Merge *all* edges that have a score below the threshold in one pass. This 
	// does not update the region graph and re-evaluate the edge scores. This is 
	// equivalent to:
	//
	//   1. find all edges below the threshold
	//   2. merge them
	//
	// This is not "associative" (in lack of a better word), i.e.:
	//
	//   merge(RAG, 4) != merge(merge(RAG, 2), 4)
	//
	// This is also too greedy. Expensive edges can be merged before cheap ones.
    for ( auto& edge: rg )
    {
		// "size" is actually the merge score
		// "weight" is probably the affinity?
        std::size_t size = func(edge.weight);

		// we don't merge 0 scores? why?
        if ( size == 0 )
        {
            break;
        }

		// get the current components of the incident regions
        ID s1 = sets.find_set(edge.id1);
        ID s2 = sets.find_set(edge.id2);

        // std::cout << s1 << " " << s2 << " " << size << "\n";

		// if the regions have not been assigned to the same component yet
		// how can s1 or s2 be zero?
		// is this a failed attempt to check if a region was involved in a merge 
		// already?
        if ( s1 != s2 && s1 && s2 )
        {
			// if the size of either region is smaller than the score computed 
			// earlier
			//
			// this implements min(size_i, size_j) < func(aff_ij)
			// for our case, func(a) = a*a*threshold
			//                     a = max_e(aff_e) for affiliated edges e 
			//                     between regions
            if ( (counts[s1] < size) || (counts[s2] < size) )
            {
				// merge the regions into s1
                counts[s1] += counts[s2];
                counts[s2]  = 0;

				// joint the connected components, the new id is s
				// (one of s1 or s2)
                ID s = sets.join(s1,s2);
				// make sure counts are correct if s == s2
                std::swap(counts[s], counts[s1]);
            }
        }
    }

	// remaps point from an old id to a new id
    std::vector<ID> remaps(counts.size());

    ID next_id = 1;

	// for each original region
    for ( ID id = 0; id < counts.size(); ++id )
    {
		// get it's connected component number
        ID s = sets.find_set(id);
		// if it has one (how could it not?)
		// and
		// has not been remapped yet
		// and
		// it is large enough (???)
        if ( s && (remaps[s] == 0) && (counts[s] >= low) )
        {
			// assign it a new id
            remaps[s] = next_id;
			// copy the counts
			// FIXME: how is this not destructive? can we guarantee that
			//   next_id < s
			// ?
			// I don't think so
            counts[next_id] = counts[s];
            ++next_id;
        }
    }

	// counts are region sizes?
	// never touched later
	// is next_id <= counts.size()?
	//   yes, we make counts smaller
    counts.resize(next_id);

    std::ptrdiff_t xdim = seg.shape()[0];
    std::ptrdiff_t ydim = seg.shape()[1];
    std::ptrdiff_t zdim = seg.shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg.data();

	// create segmentation based on remaps
    for ( std::ptrdiff_t idx = 0; idx < total; ++idx )
    {
		// TODO: this can be made faster by caching the idx
		// well, find_set is compressing the paths, so this might be okay 
		// already
        seg_raw[idx] = remaps[sets.find_set(seg_raw[idx])];
    }

    std::cout << "\tDone with remapping, total: " << (next_id-1) << std::endl;

	// create a new RG because the region ids changed and new edges have to be 
	// added
    region_graph<ID,F> new_rg;

	// a set of IDs for each region
	// are these the subsets?
    std::vector<std::set<ID>> in_rg(next_id);

	// re-visit each edge of the original region graph
    for ( auto& edge: rg )
    {
		// find the two new incident regions
        ID s1 = remaps[sets.find_set(edge.id1)];
        ID s2 = remaps[sets.find_set(edge.id2)];

		// if they are different
        if ( s1 != s2 && s1 && s2 )
        {
            auto mm = std::minmax(s1,s2);

			// if we havend seen max as neighbor of min yet
            if ( in_rg[mm.first].count(mm.second) == 0 )
            {
				// add a new edge
				// FIXME: edge weight is just the first score ever encountered! 
				// neither the max or average or whatever!
                new_rg.push_back(region_graph_edge_t<F,ID>(edge.weight, mm.first, mm.second));
                in_rg[mm.first].insert(mm.second);
            }
        }
    }

    if(recreate_rg)
        rg.swap(new_rg);

    std::cout << "\tDone with updating the region graph, size: "
              << rg.size() << std::endl;
}
