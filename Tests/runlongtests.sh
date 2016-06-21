#!/bin/bash
#set -x




cp beginreport.html testreport.html
if hash gnuplot 2>/dev/null ; then 
    GPVERSION=`gnuplot -V | awk '{print $2}'`
    echo gnuplot version $GPVERSION
    GPMAJORVERSION=`echo $GPVERSION | awk 'BEGIN {FS="."} {print $1}'`
    if [ $GPMAJORVERSION -lt 5 ] ; then
	echo You have an old gnuplot version, some test results will not display. 
    fi
else
    echo "You do not have gnuplot, many test results will not display."
    GPVERSION=0
fi


if [ ! -d Testreport ] ; then mkdir Testreport; fi
for dir in `ls`; do 
    if [ -d $dir ]; then
	if ! [ -a $dir/DONOTRUN ] ; then
	    cd $dir
	    if [ -e 'longtest.sh' ]; then
		echo running long test in $dir
		chmod +x longtest.sh
		./longtest.sh
		if [ -e 'report.html' ]; then
		    cat report.html >> ../testreport.html
		fi
	    else
		if [ -e 'run.sh' ]; then
		    echo running short test in $dir
		    chmod +x run.sh
		    ./run.sh 
		    if [ -e 'report.html' ]; then
			cat report.html >> ../testreport.html
		    fi
		fi
	    fi
	    cd ..
	fi
    fi
done

echo
echo "-----------------------------------------------------------------------------"
echo 
cat endreport.html >> testreport.html
mv testreport.html Testreport
echo The test report is in `pwd`/Testreport/testreport.html. 
echo If it does not open automatically, open it with your favorite browser. 
echo
echo "-----------------------------------------------------------------------------"
echo 


cd Testreport
if hash Open  2>/dev/null; then 
    Open testreport.html
elif hash xdg-open 2>/dev/null; then
    xdg-open testreport.html
elif hash gnome-open 2>/dev/null; then
    gnome-open testreport.html
else
    cd ..
    tar cfz tmphtmlreport.tgz Testreport
    echo "You do not appear to have a working browser".
    echo "Open the report.html file and the associated images in your favorite browser".
    echo "All the report files are compressed in tmphtmlreport.tgz." 
echo
echo "-----------------------------------------------------------------------------"
echo 
fi
