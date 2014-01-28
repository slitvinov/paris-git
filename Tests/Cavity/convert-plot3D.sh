awk '{print $2, $3, $4, $1}' ./out_$1/element-00001.dat > ./out_$1/particle-traj.3D.temp
echo ' x y z t '  > file-header.txt
cat file-header.txt ./out_$1/particle-traj.3D.temp > ./out_$1/particle-traj.3D
rm file-header.txt ./out_$1/particle-traj.3D.temp

