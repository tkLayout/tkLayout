#!/usr/bin/gawk -f

!/^<Trd1/ {
  print;
}

/^<Trd1/ {
  gsub(/^<Trd1/, "<Box");
  gsub(/ *d[xy]2="[^"]*"/, "");
  gsub(/dx1/, "dx");
  gsub(/dy1/, "dy");
  print;
}

