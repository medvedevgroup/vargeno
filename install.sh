#!/bin/bash
chmod +x ./sdsl-lite/install.sh ./sdsl-lite/build/*.sh
mkdir -p ./sdsl-lite/COMPILED
./sdsl-lite/install.sh ./sdsl-lite/COMPILED
make all
