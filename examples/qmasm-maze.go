/*
qmasm-maze generates mazes for solution by QMASM and validates solutions
returned from QMASM.

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

	"bufio"
	"errors"
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

// RandomMaze creates a maze of given dimensions.
func RandomMaze(w, h int) Maze {
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
func ColName(c int) string {
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

// RoomName maps 0-based row and column numbers to letter-number strings (e.g.,
// "J15").
func RoomName(r, c int) string {
	return fmt.Sprintf("%s%d", ColName(c), r+1)
}

// Write outputs a Maze with a given row prefix.  We assume the maze contains
// fewer than 677 columns and 1000 rows.
func (m Maze) Write(w io.Writer, pfx string) {
	// Output a column header.
	fmt.Fprint(w, pfx+"    ")
	for c := range m[0] {
		fmt.Fprintf(w, " %-2s", ColName(c))
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
# Command line:`,
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
N   0.50
E   0.50
S   0.50
W   0.50
$a1 1.00

N E   0.25
N S   0.25
N W   0.25
N $a1 0.50
E S   0.25
E W   0.25
E $a1 0.50
S W   0.25
S $a1 0.50
W $a1 0.50
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
			rstr := RoomName(r, c)
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
				fmt.Fprintf(w, "%s.N = %s.S\n", rstr, RoomName(r-1, c))
			}
			if c > 0 && !cell.W {
				fmt.Fprintf(w, "%s.W = %s.E\n", rstr, RoomName(r, c-1))
			}
		}
	}
}

// GenerateMaze generates a maze of given dimensions.
func GenerateMaze(out io.Writer, w, h int) {
	// Validate the dimensions provided.
	if w < 1 || h < 1 {
		notify.Fatal("Mazes must contain at least one row and one column")
	}
	if w >= 676 || h > 999 {
		notify.Fatal("Mazes must contain no more than 999 rows and 676 columns")
	}

	// Seed the random-number generator.
	seed := time.Now().UnixNano() * int64(os.Getpid())
	rand.Seed(seed)

	// Generate a maze.
	maze := RandomMaze(w, h)

	// Output QMASM code.
	WriteHeader(out, maze)
	WriteRooms(out, maze)
}

// NewMaze allocates a new, empty Maze.
func NewMaze() Maze {
	m := make([][]Room, 0, 10)
	for i := range m {
		m[i] = make([]Room, 0, 10)
	}
	return m
}

// Extend extends a maze by parsing a coordinate and a Boolean value.
func (m Maze) Extend(s string, b bool) Maze {
	// Parse the string.  We don't current check for validity.
	var i int
	col := 0
	for i = 0; s[i] >= 'A' && s[i] <= 'Z'; i++ {
		col = col*26 + int(s[i]-'A')
	}
	row := 0
	for ; s[i] >= '0' && s[i] <= '9'; i++ {
		row = row*10 + int(s[i]-'0')
	}
	row--
	i++ // Skip the "."
	dir := s[i]

	// Add rows and columns as necessary.  We assume we'll typically be
	// adding rooms in order.
	for row >= len(m) {
		m = append(m, make([]Room, 1))
	}
	for col >= len(m[row]) {
		m[row] = append(m[row], Room{})
	}

	// Update the specified room.
	switch dir {
	case 'N':
		m[row][col].N = b
	case 'S':
		m[row][col].S = b
	case 'E':
		m[row][col].E = b
	case 'W':
		m[row][col].W = b
	default:
		notify.Fatalf("Unknown direction \"%c\"", dir)
	}
	return m
}

// ReadSolutions returns a list of maze solutions read from a Reader.
func ReadMazes(r *bufio.Reader) []Maze {
	mazes := make([]Maze, 0, 1) // List of mazes to return
	var m Maze                  // Current maze
	for {
		// Read a line from the file and split it into fields.
		ln, err := r.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			notify.Fatal(err)
		}
		f := strings.Fields(ln)
		if len(f) == 0 {
			continue
		}

		// Start a new maze when we see "Solution".
		if f[0] == "Solution" {
			if m != nil {
				mazes = append(mazes, m)
			}
			m = NewMaze()
			continue
		}

		// Append a room to the current maze when we see a solution row.
		if len(f) == 3 && (f[2] == "True" || f[2] == "False") {
			m = m.Extend(f[0], f[2] == "True")
		}
	}
	if m != nil {
		mazes = append(mazes, m) // Final maze
	}
	return mazes
}

// NextRoom maps a "from" direction (N, S, E, or W) and room
// coordinates to a new "from" direction and room coordinates.
func (m Maze) NextRoom(r, c int, dir string) (int, int, string) {
	room := m[r][c]
	switch {
	case room.N && dir != "N":
		return r - 1, c, "S"
	case room.E && dir != "E":
		return r, c + 1, "W"
	case room.S && dir != "S":
		return r + 1, c, "N"
	case room.W && dir != "W":
		return r, c - 1, "E"
	default:
		panic("Unexpectedly stuck in a room") // Should never get here.
	}
}

// PathString returns a path through a maze or an error.
func (m Maze) PathString() (string, error) {
	// Ensure that each room has either zero or two exits.
	bool2int := map[bool]int{true: 1, false: 0}
	for r, row := range m {
		for c, cell := range row {
			n := bool2int[cell.N] + bool2int[cell.E] + bool2int[cell.S] + bool2int[cell.W]
			switch n {
			case 0:
			case 2:
			default:
				return "", fmt.Errorf("Room %s contains %d exit(s) but should have 0 or 2", RoomName(r, c), n)
			}
		}
	}

	// Find the ingress and egress columns.  Complain if more than one of
	// each is found.
	cin := -1
	for c, cell := range m[0] {
		if cell.N {
			if cin == -1 {
				cin = c
			} else {
				return "", fmt.Errorf("More than one ingress exists on the top row (%s and %d)", RoomName(0, cin), RoomName(0, c))
			}
		}
	}
	if cin == -1 {
		return "", errors.New("No ingress found in the top row")
	}
	nrows := len(m)
	cout := -1
	for c, cell := range m[nrows-1] {
		if cell.S {
			if cout == -1 {
				cout = c
			} else {
				return "", fmt.Errorf("More than one egress exists on the bottom row (%s and %d)", RoomName(nrows-1, cout), RoomName(nrows-1, c))
			}
		}
	}
	if cout == -1 {
		return "", errors.New("No egress found in the bottom row")
	}

	// Complain if the maze is open on the left or right.
	for r, row := range m {
		if row[0].W {
			return "", fmt.Errorf("Maze unexpectedly exits to the left at %s", RoomName(r, 0))
		}
		ncols := len(row)
		if row[ncols-1].E {
			return "", fmt.Errorf("Maze unexpectedly exits to the right at %s", RoomName(r, ncols-1))
		}
	}

	// Starting from the ingress, make our way to the egress.
	path := make([]string, 0, nrows*nrows)
	for r, c, d := 0, cin, "N"; r != nrows-1 || c != cout; r, c, d = m.NextRoom(r, c, d) {
		path = append(path, RoomName(r, c))
	}
	path = append(path, RoomName(nrows-1, cout))
	return strings.Join(path, " -- "), nil
}

// ValidatePaths reports if a path through a maze looks valid.
func ValidatePaths(r io.Reader) {
	// Ensure we can read line-by-line.
	rb, ok := r.(*bufio.Reader)
	if !ok {
		rb = bufio.NewReader(r)
	}

	// Read a list of mazes.
	mazes := ReadMazes(rb)

	// Process each maze in turn.
	for i, m := range mazes {
		fmt.Printf("Solution %d: ", i+1)
		s, err := m.PathString()
		if err != nil {
			fmt.Printf("[%s]\n", err)
		} else {
			fmt.Println(s)
		}
	}
}

func main() {
	// Parse the command line.
	notify = log.New(os.Stderr, path.Base(os.Args[0])+": ", 0)
	switch {
	case len(os.Args) == 4 && os.Args[1] == "gen":
		// Mode 1: Generate a maze.
		w, err := strconv.Atoi(os.Args[2])
		if err != nil {
			notify.Fatal(err)
		}
		h, err := strconv.Atoi(os.Args[3])
		if err != nil {
			notify.Fatal(err)
		}
		GenerateMaze(os.Stdout, w, h)

	case len(os.Args) == 2 && os.Args[1] == "validate":
		// Mode 2: Validate QMASM maze output.
		ValidatePaths(os.Stdin)

	default:
		notify.Fatalf("Usage: %s gen <width> <height> | %s validate", os.Args[0], os.Args[0])
	}
}
