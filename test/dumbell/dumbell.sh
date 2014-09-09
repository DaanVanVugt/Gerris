#!/bin/sh

bottom=.255
top=.49

tac <<EOF | shapes - > dumbell.gts
-0.51 -0.51
-0.51 $bottom
-0.1 $bottom
-0.1 $top
0.1 $top
0.1 $bottom
0.51 $bottom
0.51 -0.51
EOF

if gerris2D dumbell.gfs | awk '{
    if ($1 == "residual.infty:" && $3 > 6.621e-02)
      exit (1);
  }'; then
    :
else
    exit 1
fi
