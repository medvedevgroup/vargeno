#!/bin/bash
chmod +x ./sdsl-lite/install.sh ./sdsl-lite/build/*.sh
./sdsl-lite/install.sh $PREFIX
make all
