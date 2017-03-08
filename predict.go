package genexp

import (
	"bytes"
	"errors"
	"math"
	"regexp"
	"strings"
)

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

// GenerateTriplets generates the cartesian product of a string for a given
// alphabet 3 times.
func GenerateTriplets() []string {
	tripletBuffer := new(bytes.Buffer)
	reuseBuffer := new(bytes.Buffer)

	generateAllPossibleStrings(
		3, []string{"A", "T", "C", "G"}, tripletBuffer, reuseBuffer)

	tripletString := tripletBuffer.String()
	return strings.Split(
		tripletString[:len(tripletString)-len("\n")], "\n")
}

// FindAllPotentialORFs matches start "ATG" codon and stop "TGA", "TAG", "TAA" codons.
// Returns all matches, error if none.
func FindAllPotentialORFs(genome string, thresholds []int) ([]string, error) {
	stopCodons := []string{"TAA", "TAG", "TGA"}
	k := len(stopCodons[0])
	startCodonMatcher, err := regexp.Compile("ATG")
	if err != nil {
		return nil, err
	}
	stopCodonMatcher, err := regexp.Compile(strings.Join(stopCodons, "|"))
	if err != nil {
		return nil, err
	}
	matchingStartingIndices := startCodonMatcher.FindAllStringIndex(genome, -1)
	if len(matchingStartingIndices) < 1 {
		return nil, errors.New("No ORMs matched from start codons .")
	}

	matchingStopingIndices := stopCodonMatcher.FindAllStringIndex(genome, -1)
	if len(matchingStopingIndices) < 1 {
		return nil, errors.New("No ORMs matched from stop codons .")
	}

	ORFs := []string{}
	genomeLength := len(genome)
	for _, startIndex := range matchingStartingIndices {
		startOfCodon := startIndex[0]
		startOffsetByThreshold := startOfCodon + thresholds[0] + k

		// Filter out codons that are below the threshold already.
		if startOffsetByThreshold > genomeLength {
			continue
		}

		// Finding ORFs
		for _, stopIndex := range matchingStopingIndices {
			stopOfCodon := stopIndex[1]
			orfLength := stopOfCodon - startOfCodon
			if orfLength > 0 {
				if orfLength >= thresholds[0] {
					if orfLength < thresholds[1] || thresholds[1] == -1 {
						ORFs = append(ORFs, genome[startOfCodon:stopOfCodon])
					}
				}
			}
		}
	}

	if len(ORFs) < 1 {
		return nil, errors.New("No ORMs matched from start|stop codons .")
	}

	return ORFs, nil
}

// FindCodingPotential generates a codon usage table from the given ORF and
// assigns a coding potential score.
func FindCodingPotential(sequence string, codons ...string) (float64, error) {
	codonUsage, err := FindCodonUsage(sequence, codons...)
	if err != nil {
		return 0.0, err
	}
	probability := -1.0
	for codon, usage := range codonUsage {
		if usage != 0 {
			possibleUsage := math.Ceil(float64(len(sequence)) / float64(len(codon)))
			if probability == -1.0 {
				probability = (float64(usage) / possibleUsage)
			} else {
				probability *= (float64(usage) / possibleUsage)
			}
		}
	}
	if probability == -1.0 {
		return 0.0, nil
	}
	return probability, nil
}

// FindCodonUsage searches an ORF and returns a map of the number of times it
// comes up.
func FindCodonUsage(sequence string, codons ...string) (map[string]int, error) {
	codonUsage := make(map[string]int)
	for _, codon := range codons {
		codonMatcher, err := regexp.Compile(codon)
		if err != nil {
			return nil, err
		}
		codonUsage[codon] = len(codonMatcher.FindAll([]byte(sequence), -1))
	}
	return codonUsage, nil
}
