#!/bin/bash
for ((i=1;i<=10;i++))
do
  qsub ./scatter_sfr_block_linear_$i.in
done