package genexp

import (
	"bytes"
)

// insertNth inserts a newline every n-th character in a string an returns the
// new string.
func insertNth(s string, n int) string {
	buf := new(bytes.Buffer)
	prev := n - 1
	last := len(s) - 1
	for i, c := range s {
		buf.WriteRune(c)
		if i % n == prev && i != last {
			buf.WriteRune('\n')
		}
	}
	return buf.String()
}

// generateAllPossibleStrings generates all possible combinations of an alphabet
// and a given max length. Uses 2 buffers, recursively.
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

func max(n ...int) (int, int) {
  maxn, maxi := n[0], 0
  for i, a := range n {
    if a > maxn {
      maxn, maxi = a, i
    }
  }
  return maxn, maxi
}

func reverse(a string) string {
  var result string
  for i := range a {
    result = string(a[i]) + result
  }
  return result
}
