#!/bin/bash
#set -x

rm -fR out
paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

GREEN="\\033[1;32m"
NORMAL="\\033[0m"

tail -1 stats | awk '{if ($13 + $14 < 1.e-9) {print "\033[32;1m PASS\033[0m"} else {print "\033[31;1m FAIL\033[0m kinetic energy:" $13+$14}}'
