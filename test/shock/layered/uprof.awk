{
    if ($1 == "#") {
	for (i = 2; i <= NF; i++) {
	    split($i,a,":")
	    if (a[2] == "U0")
                start = a[1];
	}
    }
    else {
	dz = $5/nl;
	for (i = 0; i < nl; i++) {
	    z = dz*(0.5+i)
	    print z,$(start+i)/dz
	}
    }
}
