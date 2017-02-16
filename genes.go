package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"os"
	"path"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/olekukonko/tablewriter"
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

func main() {
	var dirPath, codons string
	var threshold int
	flag.StringVar(
		&dirPath,
		"transcripts", "transcripts", "Directory containing transcripts.")
	flag.IntVar(
		&threshold,
		"threshold", 150, "ORF Threshold, i.e. threshold = 1: ATG<A>TAG")
	flag.StringVar(
		&codons,
		"codons", "", "List of codons to determine codon usage.")
	flag.Parse()

	// Walk all files in directory and map transcript files to genome strings.
	genomes := make(map[string]string)
	filepath.Walk(dirPath, func(filePath string, info os.FileInfo, err error) error {
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
				genomes[path.Base(filePath)] = buf.String()
				return nil
			}
		}

		if err := scanner.Err(); err != nil {
			return err
		}
		return nil
	})

	var triplets []string
	if len(codons) < 1 {
		triplets = predictor.GenerateTriplets()
	} else {
		triplets = strings.Split(codons, ",")
	}

	for transcript, genome := range genomes {
		// Generate all potential ORFs from the given genome
		matches, err := predictor.FindAllPotentialORFs(genome, threshold)
		if err != nil {
			panic(err)
		}

		// Create a table to display results
		table := tablewriter.NewWriter(os.Stdout)
		table.SetHeader(
			[]string{
				"ORF",
				"CP (%)",
			},
		)

		// Find coding potential from all possible triplets for the MT genome
		// excluding the start and stop codons.
		genePredictors := make(predictor.GenePredictors, len(matches))
		for i, match := range matches {
			codingPotential, err := predictor.FindCodingPotential(
				match[3:len(match)-3], triplets...)
			if err != nil {
				panic(err)
			}
			genePredictors[i] = &predictor.GenePredictor{match, codingPotential}
		}

		// Sort by coding potential
		sort.Sort(genePredictors)

		for _, g := range genePredictors {
			// Append the information to the table for display
			table.Append(
				[]string{
					insertNth(g.Sequence, 30),
					strconv.FormatFloat(g.CodingPotential*100.0, 'E', 2, 64) + "%",
				},
			)
		}

		fmt.Println("> " + transcript)
		table.Render()
		fmt.Println()
	}
}
