# computes isoline at DRHO = 0.05
{
    while (NF == 3) {
	x = $1; 
	z = $2; r = $3;
	if (r1 > 0.05 && r < 0.05)
	    zo = z1 + (r1 - 0.05)*(z - z1)/(r1 - r);
	z1 = z; r1 = r;
	getline
    }
    print x,zo
}
