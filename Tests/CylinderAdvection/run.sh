#!/bin/bash
#set -x

rm -fR out
paris > tmpout
echo "Generating animated gif"
./gnmovie.sh
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
pariscompare out/CVoF-00000-03800.txt reference.txt $precision 1 1

hash convert 2>/dev/null || { echo >&2 "PARIS testing requires convert but it's not installed. The animated gif will not be available in
the report. Try to install Imagemagick"; exit 1; }

#
# put movie where the html report lives, if the Testreport directory exists, thus avoid creating a file called Testreport. 

if [ -d ../Testreport ] ; then
    mv CAmovie.gif ../Testreport
fi
rm -f frame*.png



GREEN="\\033[1;32m"
NORMAL="\\033[0m"
