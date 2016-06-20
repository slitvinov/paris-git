#!/bin/bash
#set -x

rm -fR out
paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
pariscompare out/CVoF-00000-03800.txt reference.txt $precision 1 1

hash convert 2>/dev/null || { echo >&2 "PARIS testing requires convert but it's not installed. The animated gif will not be available in
the report. Try to install Imagemagick"; exit 1; }

if [ $fail==0 ]; then
    echo "Generating animated gif"
    ./gnmovie.sh
fi
# put movie where the html report lives. 
mv CAmovie.gif ../Testreport
rm -f frame*.png



GREEN="\\033[1;32m"
NORMAL="\\033[0m"
