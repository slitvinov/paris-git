#! /bin/bash
set -x

# mono64-1  to mono64-512  run in less than 30 seconds on babbage cluster
# 132 Mbytes/core

make mono64-1

# 14 seconds
make mono64-8

# 21 seconds
make mono64-64

# 29 seconds
make mono64-512

make mono64-4096
make mono64-32768

