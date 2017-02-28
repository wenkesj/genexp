package align

const (
  DirectionDiagonal = iota
  DirectionHorizontal
  DirectionVertical
)

type LocalAlignment struct {
  Score int
  Alignments []string
}

type Alignment struct {
  Score int
  Direction int
}
