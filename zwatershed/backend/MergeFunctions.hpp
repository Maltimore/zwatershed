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

	MedianAffinity(const AffinityMapType& affinities) :
		_affinities(affinities) {}

	/**
	 * Get the score for an edge. An edge will be merged the earlier, the 
	 * smaller its score is.
	 */
	AffinityType operator()(const EdgeIdType& e, std::vector<EdgeIdType>& affiliatedEdges) const {

		// initial edges have their own affinity
		if (affiliatedEdges.size() == 0)
			return 1.0 - _affinities[e];

		// edges resulting from merges consult their affiliated edges

		auto median = affiliatedEdges.begin() + affiliatedEdges.size()/2;
		std::nth_element(
				affiliatedEdges.begin(),
				median,
				affiliatedEdges.end(),
				[this](const EdgeIdType& a, const EdgeIdType& b){

					return _affinities[a] < _affinities[b];
				}
		);

		return 1.0 - _affinities[*median];
	}

private:

	const AffinityMapType& _affinities;
};

#endif // MERGE_FUNCTIONS_H__

