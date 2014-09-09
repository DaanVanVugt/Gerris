if gerris2D $1 | awk '{ if ($9 < 10.) exit (1); }'; then :
    exit 1
else
    exit 0
fi
