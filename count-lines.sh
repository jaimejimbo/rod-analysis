#!/bin/bash
FILES=rods_*
for f in $FILES
do
    wc $f
done
