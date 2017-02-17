// Package predictor predicts genes from genome strings. Uses start and stop
// codons to develop open reading frames (ORFs).
package predictor

// GenePredictor is a sortable generic data structure used to identify a set of
// sequence ORFs and corresponding coding potentials.
type GenePredictor struct {
	Sequence           string
	TranscriptMatch    string
	TranscriptDistance int
	CodingPotential    float64
}

type GenePredictors []*GenePredictor

func (g GenePredictors) Len() int {
	return len(g)
}

func (g GenePredictors) Swap(i, j int) {
	g[i], g[j] = g[j], g[i]
}

func (g GenePredictors) Less(i, j int) bool {
	return g[i].TranscriptDistance > g[j].TranscriptDistance
}
