BEGIN {
  FS = "{|=|}| |[|]";
}
{
  if ($1 == "\\psfig" && $2 == "file")
    depeps = depeps " " $3;
  else if ($1 == "\\includegraphics[width" ||
	   $1 == "\\includegraphics[height" ||
	   $1 == "\\includegraphics[angle" ||
           $1 == "\\includegraphics*[width" ||
	   $1 == "\\includegraphics*[height" ||
	   $1 == "\\includegraphics*[angle") {
    for (i = 2; i <= NF; i++)
      if (match ($i, "eps"))
	depeps = depeps " " $i;
  }
}
END {
  print file "1.dvi: " depeps;
}
