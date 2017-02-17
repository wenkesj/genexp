package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"os"
	"path"
	"path/filepath"
	// "sort"
	"strconv"
	"strings"
	"sync"

	"github.com/olekukonko/tablewriter"
	"github.com/texttheater/golang-levenshtein/levenshtein"
	"github.com/wenkesj/genes/predictor"
)

// Adapted from: http://stackoverflow.com/questions/33633168/ \
// how-to-insert-a-character-every-x-characters-in-a-string-in-golang
func insertNth(s string, n int) string {
	var buf bytes.Buffer
	var prev = n - 1
	var last = len(s) - 1
	for i, c := range s {
		buf.WriteRune(c)
		if i%n == prev && i != last {
			buf.WriteRune('\n')
		}
	}
	return buf.String()
}

const (
	maxLineLength = 60
)

type Distance struct {
	distance int
	id       string
}

func minDistance(distances chan *Distance) *Distance {
	var minDistance *Distance
	for distance := range distances {
		if distance == nil {
			continue
		}
		if minDistance == nil || minDistance.distance > distance.distance {
			minDistance = distance
		}
	}
	return minDistance
}

func maxORF(array []string) string {
	maxValue := array[0]
	for _, value := range array {
		if len(value) > len(maxValue) {
			maxValue = value
		}
	}
	return maxValue
}

func main() {
	var transcriptsPath, genomePath, codons, orflengths string
	flag.StringVar(
		&transcriptsPath,
		"transcripts", "transcripts", "Directory containing transcripts.")
	flag.StringVar(
		&genomePath,
		"genome",
		"genomes/Homo_sapiens.GRCh38.dna.chromosome.MT.fa", "Path to genome.")
	flag.StringVar(
		&orflengths,
		"orflengths", "50,100",
		"ORF Threshold, i.e. orflength = 1,2: ATG<A>TAG, ATG<AA>TAG")
	flag.StringVar(
		&codons,
		"codons", "", "List of codons to determine codon usage.")
	flag.Parse()

	orflengthsString := strings.Split(orflengths, ",")
	orflengthsInt := make([]int, 2)
	if len(orflengthsString) > 1 {
		if s, err := strconv.ParseInt(orflengthsString[0], 10, 0); err == nil {
			orflengthsInt[0] = int(s)
		}
		if s, err := strconv.ParseInt(orflengthsString[1], 10, 0); err == nil {
			orflengthsInt[1] = int(s)
		}
	} else {
		if s, err := strconv.ParseInt(orflengthsString[0], 10, 0); err == nil {
			orflengthsInt[0] = int(s)
		}
		orflengthsInt[1] = -1
	}

	// Read the genome into memory
	var genome string
	file, err := os.Open(genomePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	buf := new(bytes.Buffer)
	scanner := bufio.NewScanner(file)

	// Skip the first line
	scanner.Scan()
	for scanner.Scan() {
		line := scanner.Text()
		buf.WriteString(line)
		if len(line) != maxLineLength {
			genome = buf.String()
		}
	}
	if err := scanner.Err(); err != nil {
		panic(err)
	}

	// Generate all potential ORFs from the given genome
	matches, err := predictor.FindAllPotentialORFs(genome, orflengthsInt)
	if err != nil {
		panic(err)
	}

	var triplets []string
	if len(codons) < 1 {
		triplets = predictor.GenerateTriplets()
	} else {
		triplets = strings.Split(codons, ",")
	}

	// Walk all files in directory and map transcript files to strings.
	transcripts := make(map[string]string)
	filepath.Walk(transcriptsPath, func(filePath string, info os.FileInfo, err error) error {
		if info.IsDir() {
			return nil
		}
		file, err := os.Open(filePath)
		if err != nil {
			return err
		}
		defer file.Close()

		buf := new(bytes.Buffer)
		scanner := bufio.NewScanner(file)

		// Skip the first line
		scanner.Scan()
		for scanner.Scan() {
			line := scanner.Text()
			buf.WriteString(line)
			if len(line) != maxLineLength {
				transcripts[path.Base(filePath)] = buf.String()
				return nil
			}
		}

		if err := scanner.Err(); err != nil {
			return err
		}
		return nil
	})

	// Create a table to display results
	table := tablewriter.NewWriter(os.Stdout)
	table.SetHeader(
		[]string{
			"Transcript (Distance, Lower is Better)",
			"Coding Potential (%)",
		},
	)

	// Find coding potential from all possible triplets for the MT genome
	// excluding the start and stop codons.
	fmt.Printf("Matches: %d\n", len(matches))
	genePredictors := make(predictor.GenePredictors, len(matches))
	for i, match := range matches {
		codingPotential, err := predictor.FindCodingPotential(
			match[3:len(match)-3], triplets...)
		if err != nil {
			panic(err)
		}

		// Search for transcripts and calculate the Levenshtein distance for each.
		distanceChannel := make(chan *Distance, len(transcripts))
		wait := new(sync.WaitGroup)
		for transcriptId, transcript := range transcripts {
			wait.Add(1)
			go func() {
				defer wait.Done()
				if len(transcript) > len(match) + 10 ||
						len(transcript) < len(match) - 10 {
					return
				}
				distance := levenshtein.DistanceForStrings(
					[]rune(match), []rune(transcript), levenshtein.DefaultOptions)
				distanceChannel <- &Distance{distance, transcriptId}
			}()
		}
		wait.Wait()
		close(distanceChannel)

		minTranscriptDistance := minDistance(distanceChannel)
		if minTranscriptDistance != nil {
			genePredictors = append(genePredictors,
				&predictor.GenePredictor{
					match, minTranscriptDistance.id,
					minTranscriptDistance.distance, codingPotential})
		}
		fmt.Printf("\rMatching: %.2f%%", float64(i)/float64(len(matches))*100.0)
	}
	fmt.Println()

	for _, g := range genePredictors {
		// Append the information to the table for display
		if g != nil {
			table.Append(
				[]string{
					g.TranscriptMatch + " (" +
						strconv.Itoa(g.TranscriptDistance) + ")",
					strconv.FormatFloat(g.CodingPotential*100.0, 'E', 2, 64) + "%",
				},
			)
		}
	}

	fmt.Println("> " + path.Base(genomePath))
	table.Render()
	fmt.Println()
}
