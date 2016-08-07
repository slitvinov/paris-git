#!/bin/bash
#set -x
# Perform the longtests and show the results. 

cp beginreport.html testreport.html
hash gnuplot 2>/dev/null || { echo "You do not have gnuplot, many test results will not display."; }

if [ -f Testreport ] ; then
   /bin/rm -f Testreport
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
tar cfz tmphtmlreport.tgz Testreport
echo The test report is in `pwd`/Testreport/testreport.html. 
echo All the report files are compressed in tmphtmlreport.tgz.
echo If you are testing remotely you can transfer this file to your local machine,
echo uncompress it and open testreport.html with your favorite browser.
echo 
echo "-----------------------------------------------------------------------------"

echo "Do you wish to display the report ?"
select yn in "Yes" "No"; do
    case $yn in
        Yes )     cd Testreport
    if hash Open  2>/dev/null; then 
	Open testreport.html
    elif hash xdg-open 2>/dev/null; then
	xdg-open testreport.html
    elif hash gnome-open 2>/dev/null; then
	gnome-open testreport.html
    fi
    break;;
        No ) exit;;
    esac
done
