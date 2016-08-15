# 3-input AND gate (Y = A and B and C)
!begin_macro and3
$A  -0.2727
$B   0.0000
$C  -0.2727
$Y   0.3636
$a1  0.3636

$A $B   0.0000
$A $C   0.0909
$A $Y  -0.1818
$A $a1 -0.1818
$B $C   0.0000
$B $Y  -0.3636
$B $a1  0.3636
$C $Y  -0.1818
$C $a1 -0.1818
$Y $a1 -0.1818
!end_macro and3

# 2-input OR gate (Y = A or B)
!begin_macro or2
$A  0.3333
$B  0.3333
$Y -0.6667

$A $B  0.3333
$A $Y -0.6667
$B $Y -0.6667
!end_macro or2

# 1-input NOT gate (Y = not A)
!begin_macro not1
$A -0.5
$Y -0.5

$A $Y 1.0
!end_macro not1
