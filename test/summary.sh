endstamp=`date +%s`
end=`date -d@$endstamp +"%a %d %b %H:%M:%S"`

dirs="$TESTS"

status=`for d in $dirs; do find $d -name status; done`

fail=`grep FAIL $status | wc -l`
pass=`grep PASS $status | wc -l`
n=`echo $status | wc -w`

logs=""
for d in $dirs; do
    logs=$d.sh.log
done

earliest=`ls -tr *.log | head -n1`
startstamp=`stat -c %Y $earliest`
start=`date -d@$startstamp +"%a %d %b %H:%M:%S"`

system=`uname -o -n -m | sed 's/_/\\\\_/'`
path=`which gerris2D`
version=`gerris2D -V 2>&1 | head -1 | cut -d' ' -f6-`

if test x$fail = x0; then
    status="{\color{OliveGreen}PASS ($pass)}"
else
    status="{\color{Red}FAIL ($fail/$n)}"
fi

elapsed=`awk "BEGIN{
  elapsed = $endstamp - $startstamp
  days = int(elapsed/86400)
  elapsed -= days*86400
  hours = int(elapsed/3600)
  elapsed -= hours*3600
  mins = int(elapsed/60)
  elapsed -= mins*60
  if (days > 0)
    printf(\"%02d:\", days);
  if (hours > 0)
    printf(\"%02d:\", hours);
  printf(\"%02d:%02d\", mins, elapsed);
}"`

cat <<EOF > summary.tex
\begin{tabular}{ll}
{\bf Version} & $version \\\\
{\bf Path} & $path \\\\
{\bf System} & $system \\\\
{\bf Start} & $start \\\\
{\bf Finish} & $end \\\\
{\bf Elapsed} & $elapsed \\\\
{\bf Status} & $status
\end{tabular}
EOF
