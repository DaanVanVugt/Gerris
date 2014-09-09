#!/bin/sh

for module in `cat modules.list`; do
    for dim in 2D 3D; do
	lib=$1/.libs/lib"$module""$dim".so
	if test ! -f $lib; then
	    lib=$2/lib"$module""$dim".so
	fi
	if test -f $lib; then
	    nm -fb $lib | grep ".* T gfs_.*_class$" | grep -v "gfs_gl" | cut -d" " -f3-4
	fi
    done
done | sort | uniq | sed -e 's/_class//g' -e 's/^./\U&/' -e 's/_./\U&/g' -e 's/_//g' | \
awk '{ print $0; gsub ("^Gfs", ""); print $0; }'
