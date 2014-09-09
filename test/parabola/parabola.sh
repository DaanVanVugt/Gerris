levels="5 6 7 8 9"

if test x$donotrun != xtrue; then
    for level in $levels; do
	if gerris2D -DLEVEL=$level $1; then :
	else
	    exit 1
	fi
    done
fi

for level in $levels; do
    awk -v level=$level 'BEGIN{
      s1 = 0.;
      s2 = 0.;
      smax = 0.;
      n = 0;
      h0 = 10.;
    }{
      n++;
      s1 += $5;
      s2 += $7*$7;
      if ($9 > smax) smax = $9;
    }END { print level, s1/n/h0, sqrt(s2/n)/h0, smax/h0; }' < error-$level
done > error

for level in $levels; do
    if  paste U-$level vol-$level | awk -v level=$level 'BEGIN {
      sum = 0.; 
      n = 0; 

      h0 = 10.
      a = 3000.
      tau = 1e-3
      B = 5.
      G = 9.81
      p = sqrt (8.*G*h0)/a
      s = sqrt (p*p - tau*tau)/2.
    } {
          u0 = $5/$10;
          t = $3;
          ref = B*exp (-tau*t/2.)*sin (s*t);
          sum += (u0 - ref)*(u0 - ref);
          n += 1;
        }
        END {
          print level, sqrt (sum/n);
        }'; then
	:
    else
	exit 1;
    fi
done > error-u

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20 enhanced

    h0 = 10.
    a = 3000.
    tau = 1e-3
    B = 5.
    G = 9.81
    p = sqrt (8.*G*h0)/a
    s = sqrt (p*p - tau*tau)/2.
    u0(t) = B*exp (-tau*t/2.)*sin (s*t)

    set output 'convergence.eps'
    set xlabel 'Level'
    set ylabel 'Relative error norms'
    set key bottom left
    set logscale y
    set xtics 0,1
    set grid
    ftitle(a,b) = sprintf("order %4.2f", -b/log(2.))
    f1(x)=a1+b1*x
    fit f1(x) 'error' u 1:(log(\$2)) via a1,b1
    f2(x)=a2+b2*x
    fit f2(x) 'error' u 1:(log(\$3)) via a2,b2
    fm(x)=am+bm*x
    fit fm(x) 'error' u 1:(log(\$4)) via am,bm
    plot 'error.ref' u 1:2 t '|h|_1 (ref)' ps 1.5, \
         'error.ref' u 1:3 t '|h|_2 (ref)' ps 1.5, \
         'error.ref' u 1:4 t '|h|_{max} (ref)' ps 1.5, \
         exp (f1(x)) t ftitle(a1,b1), \
         exp (f2(x)) t ftitle(a2,b2), \
         exp (fm(x)) t ftitle(am,bm),  \
         'error' u 1:2 t '|h|_1' ps 1.5, \
         'error' u 1:3 t '|h|_2' ps 1.5, \
         'error' u 1:4 t '|h|_{max}' ps 1.5

    set output 'convergence-u.eps'
    set xlabel 'Level'
    set ylabel '|u_0|_2'
    set key top right
    set logscale y
    set xtics 0,1
    set grid
    fit f2(x) 'error-u' u 1:(log(\$2/B)) via a2,b2
    plot exp (f2(x)) t ftitle(a2,b2), \
         'error-u' u 1:(\$2/B) t '' ps 1.5

    set output 'u0.eps'

    set xtics auto
    set ytics auto
    unset grid
    unset logscale
    set key top right
    set ylabel 'u0'
    set xlabel 'Time'
    plot u0(x) t 'Analytical', '< paste U-7 vol-7' u 3:(\$5/\$10) every 2 w p t 'Numerical'

    set output 'elevation.eps'
    set xlabel 'x (m)'
    set ylabel 'z (m)'
    t = 1500
    psi(x) = a*a*B*B*exp (-tau*t)/(8.*G*G*h0)*(- s*tau*sin (2.*s*t) + \
      (tau*tau/4. - s*s)*cos (2.*s*t)) - B*B*exp(-tau*t)/(4.*G) - \
      exp (-tau*t/2.)/G*(B*s*cos (s*t) + tau*B/2.*sin (s*t))*x + h0
    bed(x) = h0*(x/a)**2
    set key top center
    plot [-5000:5000] \
      'sim-6-1500.txt' u 1:7:8 w filledcu lc 3 t 'Numerical', \
      psi(x) > bed(x) ? psi(x) : bed(x) lc 2 t 'Analytical', \
      bed(x) lw 3 lc 1 lt 1 t 'Bed profile'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error',1,4) - Curve('error.ref',1,4)).max() > 1e-4:
    print (Curve('error',1,4) - Curve('error.ref',1,4)).max()
    exit(1)
EOF
else
   exit 1
fi
