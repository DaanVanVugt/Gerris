function replace_params(s, b,    i)
{
    for (i in b)
	gsub(b[i], "($" i ")", s);
    return s;
}

BEGIN {
    print prefix "changecom()" prefix "dnl";
}
{
    if ($1 == "GfsDefine" || $1 == "Define") {
	macro = $2;
	split("", b); # delete b

	# if we could use gawk, we could do ...
#	if (match(macro, /(.+)\((.+)\)/, a)) {
#	    macro = a[1];
#	    split(a[2],b,",");
#	}
	# but for portability we need to do ...
	if (match(macro, /.+\(/)) {
	    a[1] = substr(macro, RSTART, RLENGTH - 1);
	    last = substr(macro, RSTART + RLENGTH);
	    if (match(last, /.+\)/)) {
		a[2] = substr(last, RSTART, RLENGTH - 1);
		macro = a[1];
		split(a[2],b,",");
	    }
	}
	printf (prefix "define(`%s',`%s", macro, replace_params($3, b));
	for (i = 4; i <= NF; i++)
	    printf (" %s", replace_params($i, b));
	printf ("')\n");
    }
    else if ($1 == "GfsInclude" || $1 == "Include")
	printf (prefix "include(%s)\n", $2);
    else
	print $0;
}
