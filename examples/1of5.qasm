################################################
# QASM example: one "on" bit out of five total #
# By Scott Pakin <pakin@lanl.gov>              #
################################################

# Point weights
# A-E are the variables we care about.
# $a1-$a3 are ancillary variables.
A   -2
B    1
C   -2
D   -1
E    1
$a1  0
$a2  4
$a3  3

# Coupler strengths
A   B    2
A   C    2
A   D    2
A   E    1
A   $a1  1
A   $a2 -4
A   $a3 -4
B   C    2
B   D    2
B   E    1
B   $a1  1
B   $a2 -4
B   $a3 -1
C   D    2
C   E    3
C   $a1 -1
C   $a2 -4
C   $a3 -4
D   E    4
D   $a1 -2
D   $a2 -4
D   $a3 -3
E   $a1 -4
E   $a2  0
E   $a3 -4
$a1 $a2 -4
$a1 $a3  3
$a2 $a3  0
