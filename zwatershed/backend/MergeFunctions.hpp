#ifndef MERGE_FUNCTIONS_H__
#define MERGE_FUNCTIONS_H__

/**
 * Scores edges with 1 - median affinity.
 */
template <typename NodeIdType, typename AffinityType>
class MedianAffinity {

public:

	typedef RegionGraph<NodeIdType>                                  RegionGraphType;
	typedef typename RegionGraphType::template EdgeMap<AffinityType> AffinityMapType;
	typedef typename RegionGraphType::EdgeIdType                     EdgeIdType;

	MedianAffinity(RegionGraphType& regionGraph, const AffinityMapType& affinities) :
		_affinities(affinities),
		_affiliatedEdges(regionGraph) {}

	/**
	 * Get the score for an edge. An edge will be merged the earlier, the 
	 * smaller its score is.
	 */
	AffinityType operator()(EdgeIdType e) {

		std::vector<EdgeIdType>& affiliatedEdges = _affiliatedEdges[e];

		// initial edges have their own affinity
		if (affiliatedEdges.size() == 0)
			return 1.0 - _affinities[e];

		// edges resulting from merges consult their affiliated edges

		auto median = affiliatedEdges.begin() + affiliatedEdges.size()/2;
		std::nth_element(
				affiliatedEdges.begin(),
				median,
				affiliatedEdges.end(),
				[this](EdgeIdType a, EdgeIdType b){

					return _affinities[a] < _affinities[b];
				}
		);

		return 1.0 - _affinities[*median];
	}

	void notifyEdgeMerge(EdgeIdType from, EdgeIdType to) {

		if (_affiliatedEdges[from].size() == 0)
			// 'from' is an initial edge
			_affiliatedEdges[to].push_back(from);
		else
			// 'from' is a compound edge, copy its affiliated edges to the new 
			// affiliated edge list
			std::copy(
					_affiliatedEdges[from].begin(),
					_affiliatedEdges[from].end(),
					std::back_inserter(_affiliatedEdges[to]));

		// clear affiliated edges of merged region edges -- they are not 
		// needed anymore
		_affiliatedEdges[from].clear();
	}

private:

	const AffinityMapType& _affinities;

	// for every new edge between regions u and v, the edges of the initial RAG 
	// between any child of u and any child of v
	//
	// initial edges will have this empty
	typename RegionGraphType::template EdgeMap<std::vector<EdgeIdType>> _affiliatedEdges;
};

#endif // MERGE_FUNCTIONS_H__

