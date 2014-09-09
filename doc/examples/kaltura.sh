# insert the Kaltura.org javascript for HTML5 video
# http://www.kaltura.org/project/HTML5_Video_Media_JavaScript_Library

f=`echo $1 | sed 's/\.html$//'`
if grep -q "</video>" $f.html; then
    tmp=`mktemp /tmp/kaltura.XXXXXXXXXX`
    sed 's/<\/HEAD>/<script type="text\/javascript" src="http:\/\/html5.kaltura.org\/js"><\/script>\n<\/HEAD>/' < $f.html > $tmp
    mv -f $tmp $f.html
fi
