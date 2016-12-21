/*
generate-maze generates a maze for solution by QMASM.

Author: Scott Pakin, pakin@lanl.gov
*/
package main

import (
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"path"
	"strconv"
	"strings"
	"time"

	"github.com/spakin/disjoint"
)

// notify is used to write status messages to the user.
var notify *log.Logger

// A Room is identified by its walls and by the other rooms it can reach.
type Room struct {
	N       bool              // North side of room is a wall
	S       bool              // South side of room is a wall
	E       bool              // East side of room is a wall
	W       bool              // West side of room is a wall
	Reaches *disjoint.Element // Element in a set of reachable rooms
}

// A Maze is a 2-D array of Rooms.
type Maze [][]Room

// NewMaze creates a maze of given dimensions.
func NewMaze(w, h int) Maze {
	// Allocate and initialize the maze to all walls present.
	maze := make([][]Room, h)
	for y := range maze {
		maze[y] = make([]Room, w)
		for x := range maze[y] {
			// Start with all walls present and no other rooms reachable.
			maze[y][x].N = true
			maze[y][x].S = true
			maze[y][x].E = true
			maze[y][x].W = true
			maze[y][x].Reaches = disjoint.NewElement()
		}
	}

	// Repeatedly remove walls until a single connected component remains.
	for cc := w * h; cc > 1; {
		// Because of symmetry, we need only connect to the right or
		// down rather than in all four directions.
		x0 := rand.Intn(w)
		y0 := rand.Intn(h)
		x1 := x0
		y1 := y0
		dir := rand.Intn(2)
		if dir == 0 && x0 < w-1 {
			x1++ // Go right.
		} else if dir == 1 && y0 < h-1 {
			y1++ // Go down.
		} else {
			continue // Can't go in the desired direction
		}
		if maze[y0][x0].Reaches.Find() == maze[y1][x1].Reaches.Find() {
			continue // Already connected
		}

		// Tear down the wall.
		if dir == 0 {
			// Right/left
			maze[y0][x0].E = false
			maze[y1][x1].W = false
		} else {
			// Down/up
			maze[y0][x0].S = false
			maze[y1][x1].N = false
		}
		disjoint.Union(maze[y0][x0].Reaches, maze[y1][x1].Reaches)
		cc--
	}

	// Punch one hole on the top and one hole on the bottom.
	maze[0][rand.Intn(w)].N = false
	maze[h-1][rand.Intn(w)].S = false

	// Return the generated maze.
	return maze
}

// ColName maps the numbers {0, 1, 2, ..., 675} to the strings {"A", "B", "C",
// ..., "ZZ"}.
func (m Maze) ColName(c int) string {
	switch {
	case c < 26:
		return fmt.Sprintf("%c", c+'A')

	case c < 26*26:
		return fmt.Sprintf("%c%c", c/26+'A'-1, c%26+'A')

	default:
		notify.Fatalf("Invalid column number %d", c)
	}
	return "" // Will never get here.
}

// Write outputs a Maze with a given row prefix.  We assume the maze contains
// fewer than 677 columns and 1000 rows.
func (m Maze) Write(w io.Writer, pfx string) {
	// Output a column header.
	fmt.Fprint(w, pfx+"    ")
	for c := range m[0] {
		fmt.Fprintf(w, " %-2s", m.ColName(c))
	}
	fmt.Fprint(w, "\n")

	// Output each row in turn.
	for r, row := range m {
		// Output the row's northern walls as one line of ASCII graphics.
		fmt.Fprint(w, pfx+"    ")
		for _, cell := range row {
			if cell.N {
				fmt.Fprintf(w, "+--")
			} else {
				fmt.Fprintf(w, "+  ")
			}
		}
		fmt.Fprintln(w, "+")

		// Output the row's western walls as another line of ASCII
		// graphics.
		fmt.Fprintf(w, "%s%3d ", pfx, r+1)
		for _, cell := range row {
			if cell.W {
				fmt.Fprintf(w, "|  ")
			} else {
				fmt.Fprintf(w, "   ")
			}
		}

		// End the line with the single, easternmost wall.
		if row[len(row)-1].E {
			fmt.Fprintln(w, "|")
		} else {
			fmt.Fprintln(w, "")
		}
	}

	// Output the bottomost row's southern walls as a final line of ASCII
	// graphics.
	fmt.Fprint(w, pfx+"    ")
	for _, cell := range m[len(m)-1] {
		if cell.S {
			fmt.Fprintf(w, "+--")
		} else {
			fmt.Fprintf(w, "+  ")
		}
	}
	fmt.Fprintln(w, "+")
}

// WriteHeader writes some boilerplate header to a file.
func WriteHeader(w io.Writer, m Maze) {
	fmt.Fprintln(w, `#########################################
# Find the shortest path through a maze #
# By Scott Pakin <pakin@lanl.gov>       #
#########################################

# This is a generated file.
# Command line: `,
		strings.Join(os.Args, " "), `

# Maze to solve:
#`)
	m.Write(w, "# ")
	fmt.Fprintln(w, `
# Truth table for a room:
#
#   0 0 0 0
#   0 0 1 1
#   0 1 0 1
#   0 1 1 0
#   1 0 0 1
#   1 0 1 0
#   1 1 0 0

# Define a macro for a room that has the preceding truth table as
# the degenerate ground state of the corresponding Hamiltonian.
!begin_macro room
N    0.25
E    0.00
S    0.00
W    0.25
$a1  1.00
$a2  1.00

N   E    0.50
N   S    0.50
N   W    0.50
N   $a1  1.00
N   $a2 -0.75
E   S    0.50
E   W    0.50
E   $a1  1.00
E   $a2 -1.00
S   W    0.50
S   $a1  1.00
S   $a2 -1.00
W   $a1  1.00
W   $a2 -0.75
$a1 $a2 -1.00
!end_macro room

# Define some helpful aliases.
!alias egress TRUE
!alias wall   FALSE
`)
}

// WriteRooms writes the rooms of a maze to a file.
func WriteRooms(w io.Writer, m Maze) {
	fmt.Fprintln(w, "# Output in turn each room of the maze.")
	for r, row := range m {
		for c, cell := range row {
			fmt.Fprintln(w, "")
			rstr := fmt.Sprintf("%s%d", m.ColName(c), r+1)
			fmt.Fprintf(w, "!use_macro room %s\n", rstr)

			// Output egresses (which always face north or south in
			// the current implementation).
			if r == 0 && !cell.N {
				fmt.Fprintln(w, rstr+".N := egress")
			}
			if r == len(m)-1 && !cell.S {
				fmt.Fprintln(w, rstr+".S := egress")
			}

			// Output all walls.
			if cell.N {
				fmt.Fprintln(w, rstr+".N := wall")
			}
			if cell.E {
				fmt.Fprintln(w, rstr+".E := wall")
			}
			if cell.S {
				fmt.Fprintln(w, rstr+".S := wall")
			}
			if cell.W {
				fmt.Fprintln(w, rstr+".W := wall")
			}

			// Output northern and western paths.  (The others are
			// symmetric.)
			if r > 0 && !cell.N {
				fmt.Fprintf(w, "%s.N = %s%d.S\n", rstr, m.ColName(c), r)
			}
			if c > 0 && !cell.W {
				fmt.Fprintf(w, "%s.W = %s%d.E\n", rstr, m.ColName(c-1), r+1)
			}
		}
	}
}

func main() {
	// Parse the command line.
	notify = log.New(os.Stderr, path.Base(os.Args[0])+": ", 0)
	if len(os.Args) != 3 {
		notify.Fatalf("Usage: %s <width> <height>", os.Args[0])
	}
	width, err := strconv.Atoi(os.Args[1])
	if err != nil {
		notify.Fatal(err)
	}
	height, err := strconv.Atoi(os.Args[2])
	if err != nil {
		notify.Fatal(err)
	}
	if width < 1 || height < 1 {
		notify.Fatal("Mazes must contain at least one row and one column")
	}
	if width >= 676 || height > 999 {
		notify.Fatal("Mazes must contain no more than 999 rows and 676 columns")
	}

	// Seed the random-number generator.
	seed := time.Now().UnixNano() * int64(os.Getpid())
	rand.Seed(seed)

	// Generate a maze.
	maze := NewMaze(width, height)

	// Output QMASM code.
	WriteHeader(os.Stdout, maze)
	WriteRooms(os.Stdout, maze)
}
