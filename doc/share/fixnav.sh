for f in $1/*.html; do
    sed 's/contents_motif.gif/contents.png/g' < $f | \
    sed 's/next_motif.gif/next.png/g' | \
    sed 's/previous_motif.gif/prev.png/g' \
    > $f.bak
    mv -f $f.bak $f
done
