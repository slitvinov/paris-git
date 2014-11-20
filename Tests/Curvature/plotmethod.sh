ndepth=`head -60  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
sed s/ndepth/$ndepth/g plotmethod.gp > pm.tmp
gnuplot < pm.tmp
