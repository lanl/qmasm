# Solve a circuit-satisfiability problem.

!include <gates>

!use_macro not1 not_x4
not_x4.$A = x3
not_x4.$Y = $x4

!use_macro or2 or_x5
or_x5.$A = x1
or_x5.$B = x2
or_x5.$Y = $x5

!use_macro not1 not_x6
not_x6.$A = $x4
not_x6.$Y = $x6

!use_macro and3 and_x7
and_x7.$A = x1
and_x7.$B = x2
and_x7.$C = $x4
and_x7.$Y = $x7

!use_macro or2 or_x8
or_x8.$A = $x5
or_x8.$B = $x6
or_x8.$Y = $x8

!use_macro or2 or_x9
or_x9.$A = $x6
or_x9.$B = $x7
or_x9.$Y = $x9

!use_macro and3 and_x10
and_x10.$A = $x8
and_x10.$B = $x9
and_x10.$C = $x7
and_x10.$Y = x10
