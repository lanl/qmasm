####################################
# Comparator for a sorting network #
# By Scott Pakin <pakin@lanl.gov>  #
####################################

# Semantics:
#
#   IF $a < $b THEN
#     $min = $a
#     $max = $b
#   ELSE
#     $min = $b
#     $max = $a
#   ENDIF

!begin_macro comparator
$a    0
$b    0 
$min  1
$max -1

$a $b      1
$a $min   -1
$a $max   -0.5
$b $min   -1
$b $max   -0.5
$min $max -0.5
!end_macro comparator
