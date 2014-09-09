for f in height*.gfs; do
    if gerris2D $f > ref.gfs; then :
    else
	echo "  FAIL: gerris2D $f"
	exit 1
    fi

    np=`awk 'BEGIN{max = 0}{ 
               if ($3 == "pid" && $4 == "=" && $5 > max) 
                 max = $5;
             }END{print max+1}' < $f`
    if mpirun -np $np gerris2D $f > run.gfs; then :
    else
	echo "  FAIL: mpirun -np $np gerris2D $f"
	exit 1
    fi
    
    for v in T_Hbx T_Hby T_Htx T_Hty K; do
	if gfscompare2D -v run.gfs ref.gfs $v 2> log; then :
	else
	    cat log
	    echo "  FAIL: $v $f"
	    exit 1
	fi
	if awk '{ if ($1 == "total" && $8 > 1e-10) exit 1; }' < log; then :
	else
	    cat log
	    echo "  FAIL: $v $f"
	    exit 1
	fi
    done
done
