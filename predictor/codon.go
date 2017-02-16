package predictor

import (
	"bytes"
	"errors"
	"math"
	"regexp"
	"strings"

	"github.com/fatih/color"
)

func generateAllPossibleStrings(
	maxLength int, alphabet []string,
	tripletBuffer, reuseBuffer *bytes.Buffer) {
	if reuseBuffer.Len() == maxLength {
		tripletBuffer.WriteString(reuseBuffer.String() + "\n")
	} else {
		for _, char := range alphabet {
			oldCurr := bytes.NewBufferString(reuseBuffer.String())
			reuseBuffer.WriteString(char)
			generateAllPossibleStrings(maxLength, alphabet, tripletBuffer, reuseBuffer)
			reuseBuffer.Reset()
			reuseBuffer.WriteString(oldCurr.String())
		}
	}
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

// maxIndex finds the greatest valued integer from a 2D slice of index ranges.
func maxIndex(array [][]int) int {
	maxValue := array[0][0]
	for _, subArray := range array {
		value := subArray[0]
		if value > maxValue {
			maxValue = value
		}
	}
	return maxValue
}

func maxORFs(array []string) []string {
	maxs := []string{}
	maxValue := array[0]
	for _, value := range array {
		if len(value) > len(maxValue) {
			maxValue = value
		}
	}

	for _, value := range array {
		if len(value) == len(maxValue) {
			maxs = append(maxs, value)
		}
	}
	return maxs
}

// FindAllPotentialORFs matches start "ATG" codon and stop "TGA", "TAG", "TAA" codons.
// Returns all matches, error if none.
func FindAllPotentialORFs(genome string, threshold int) ([]string, error) {
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
		return nil, errors.New("No ORMs matched from start|stop codons .")
	}

	ORFs := []string{}
	genomeLength := len(genome)
	for _, startIndex := range matchingStartingIndices {
		startOfCodon := startIndex[1]
		startOffsetByThreshold := startOfCodon + threshold + k
		// Check if the next possible ORF is at the end, if so quit
		if startOffsetByThreshold-genomeLength > 0 {
			startOffsetByThreshold = genomeLength - 1
		}

		// Starting ORF from sequence
		sequence := genome[startOfCodon:startOffsetByThreshold]
		matchingStopingIndices := stopCodonMatcher.FindAllStringIndex(sequence, -1)
		if len(matchingStopingIndices) < 1 {
			continue
		}
		ORFs = append(ORFs, sequence[:maxIndex(matchingStopingIndices)+k])
	}

	if len(ORFs) < 1 {
		return nil, errors.New("No ORMs matched from start|stop codons .")
	}

	return maxORFs(ORFs), nil
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

// ColorSequence colors a target string within a sequence a certain color.
// Returns the colored string.
func ColorSequence(c, sequence string, targets ...string) string {
	var colorString func(a ...interface{}) string
	switch strings.ToLower(c) {
	case "blue":
		colorString = color.New(color.FgBlue).SprintFunc()
		break
	case "yellow":
		colorString = color.New(color.FgYellow).SprintFunc()
		break
	case "red":
		colorString = color.New(color.FgRed).SprintFunc()
		break
	case "green":
		colorString = color.New(color.FgGreen).SprintFunc()
		break
	default:
		colorString = color.New(color.FgBlue).SprintFunc()
	}
	var copySequence string
	for _, target := range targets {
		copySequence = strings.Replace(sequence, target, colorString(target), -1)
	}
	return copySequence
}
