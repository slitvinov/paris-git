#!/bin/bash
#set -x

rm -fR out
paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

compare out/CVoF-00000-03800.txt reference.txt $precision 1 1

GREEN="\\033[1;32m"
NORMAL="\\033[0m"



