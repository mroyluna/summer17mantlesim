#!/bin/bash
sed -e '/\<Points\>/,/\<\/Points\>/!d' $1 | perl -pe 's/\<.*?\>//g;' -e 's/\ \ /\n/g;' > Points
sed -e '/\<PointData/,/\<\/PointData/!d' $1 | perl -pe 's/\<.*?\>//g;' -e 's/\ \ /\n/g;' > Data
pr -mt -s\    Points Data | awk '{print $1, $2, $4, $5}' > PythonSoln
rm Points
rm Data
