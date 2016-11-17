#ifndef MERGE_FUNCTIONS_H__
#define MERGE_FUNCTIONS_H__

/**
 * Scores edges with 1 - median affinity.
 */
class MedianAffinity {

public:

	/**
	 * Get the score for an edge. An edge will be merged the earlier, the 
	 * smaller its score is.
	 */
	template <typename EdgeType>
	typename EdgeType::AffinityType operator()(const EdgeType& edge, std::vector<EdgeType>& affiliatedEdges) const {

		// initial edges have their own affinity
		if (affiliatedEdges.size() == 0)
			return 1.0 - edge.affinity;

		// edges resulting from merges consult their affiliated edges

		auto median = affiliatedEdges.begin() + affiliatedEdges.size()/2;
		std::nth_element(
				affiliatedEdges.begin(),
				median,
				affiliatedEdges.end(),
				[](const EdgeType& a, const EdgeType& b){

					return a.affinity < b.affinity;
				}
		);

		return 1.0 - median->affinity;
	}
};

#endif // MERGE_FUNCTIONS_H__

