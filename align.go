package genexp

import (
	"fmt"
	"bytes"
)

type alignmentDirection int

const (
	vertical alignmentDirection = iota
	diagonal
	horizontal
	none
)

type tableAlignment struct {
	Score     int
	Direction alignmentDirection
}

// LocalAlignmentCost is the cost function for local alignments.
func LocalAlignmentCost(a, b rune) int {
	if a == b {
		return 1
	}
	return -1
}

// Align computes the optimal local alignment of two sequences. Sequence a and
// b are strings. Cost is a function to calculate the cost based on given type.
//
// Align returns all optimal local alignments for a and b and their
// corresponding score.
func Align(a, b string, cost func(rune, rune) int) ([]string, []string, int) {
	var align func(*tableAlignment, *bytes.Buffer, *bytes.Buffer, int, int)
	rows := len(a) + 1
	columns := len(b) + 1
	dTable := make([][]*tableAlignment, rows)
	for i := range dTable {
		dTable[i] = make([]*tableAlignment, columns)
	}

	dTable[0][0] = &tableAlignment{
		Score:     0,
		Direction: alignmentDirection(none),
	}

	for i := range dTable[1:] {
		dTable[i+1][0] = &tableAlignment{
			Score:     cost(rune(a[i]), 0),
			Direction: vertical,
		}
	}
	for j := range dTable[0][1:] {
		dTable[0][j+1] = &tableAlignment{
			Score:     cost(0, rune(b[j])),
			Direction: horizontal,
		}
	}

	fmt.Print("    ")
	for _, c := range b {
		fmt.Print(string(c) + " ")
	}
	fmt.Println()

	fmt.Print("  0 ")
	i := 0
	for i < len(b) {
		fmt.Print("0 ")
		i++
	}
	fmt.Println()

	maxAlignments := []*tableAlignment{dTable[0][0]}
	maxRows := []int{0}
	maxColumns := []int{0}
	for row := range dTable[:rows-1] {
		fmt.Print(string(a[row]) + " 0 ")
		for column := range dTable[row][:columns-1] {
			nextRow := row + 1
			nextColumn := column + 1
			score, direction := max(
				dTable[nextRow][column].Score+cost(rune(a[row]), 0),
				dTable[row][column].Score+cost(rune(a[row]), rune(b[column])),
				dTable[row][nextColumn].Score+cost(0, rune(b[column])),
				0,
			)
			alignment := &tableAlignment{
				Score:     score,
				Direction: alignmentDirection(direction),
			}
			dTable[nextRow][nextColumn] = alignment

			if alignment.Score > maxAlignments[0].Score {
				maxAlignments = []*tableAlignment{alignment}
				maxRows = []int{nextRow}
				maxColumns = []int{nextColumn}
			} else if alignment.Score == maxAlignments[0].Score {
				maxAlignments = append(maxAlignments, alignment)
				maxRows = append(maxRows, nextRow)
				maxColumns = append(maxColumns, nextColumn)
			}
			fmt.Print(alignment.Score, " ")
		}
		fmt.Println()
	}

	align = func(alignment *tableAlignment, aBuf, bBuf *bytes.Buffer, i, j int) {
		if alignment.Direction == none {
			if i-1 >= 0 && j-1 >= 0 {
				align(dTable[i-1][j-1], aBuf, bBuf, i-1, j-1)
			} else {
				aBuf.WriteString("...")
				bBuf.WriteString("...")
			}
		} else if alignment.Direction == vertical {
			aBuf.WriteString(string(a[i-1]))
			bBuf.WriteString("_")
			align(dTable[i-1][j], aBuf, bBuf, i-1, j)
		} else if alignment.Direction == horizontal {
			aBuf.WriteString("_")
			bBuf.WriteString(string(b[j-1]))
			align(dTable[i][j-1], aBuf, bBuf, i, j-1)
		} else if alignment.Direction == diagonal {
			aBuf.WriteString(string(a[i-1]))
			bBuf.WriteString(string(b[j-1]))
			align(dTable[i-1][j-1], aBuf, bBuf, i-1, j-1)
		}
	}

	aAlignments, bAlignments := make([]string, len(maxAlignments)), make([]string, len(maxAlignments))
	for i, maxAlignment := range maxAlignments {
		aBuf, bBuf := new(bytes.Buffer), new(bytes.Buffer)
		align(maxAlignment, aBuf, bBuf, maxRows[i], maxColumns[i])
		aRev := reverse(aBuf.String())
		bRev := reverse(bBuf.String())
		aAlignments[i], bAlignments[i] = aRev+"...", bRev+"..."
	}
	return aAlignments, bAlignments, maxAlignments[0].Score
}
