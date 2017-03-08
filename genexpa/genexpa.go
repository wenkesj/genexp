package main

import (
	"fmt"
	"os"
	"github.com/wenkesj/genexp"
)

func main() {
	if len(os.Args) < 3 {
		panic("Expected 2 sequences as part of the command. \n For example: genexpa ATG GTG")
	}
	a := os.Args[1]
	b := os.Args[2]

	alignedAs, alignedBs, score := genexp.Align(a, b, genexp.LocalAlignmentCost)
	fmt.Println("Alignments: ")
	for i := range alignedAs {
		fmt.Println("A: ", alignedAs[i])
		fmt.Println("B: ", alignedBs[i])
		fmt.Println("Score: ", score)
	}
}
