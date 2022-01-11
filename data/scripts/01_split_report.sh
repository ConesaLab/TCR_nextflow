#!/bin/bash
for file in ${params.outdir}/*.report; do
   csplit -f $file $file '/^==/' '{*}'
   rm *.report05
done
