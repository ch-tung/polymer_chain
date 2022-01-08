#!/bin/bash
for ((i=1;i<=9;i++))
do
  cp ./scatter_sfr_block_linear_0.py ./scatter_sfr_block_linear_$i.py
  sed -i "40s/i_s0 = 0/i_s0 = $i/g" ./scatter_sfr_block_linear_$i.py
  
  cp ./scatter_sfr_block_linear_0.in ./scatter_sfr_block_linear_$i.in
  sed -i "3s/scatter_sfr_block_linear_0/scatter_sfr_block_linear_$i/g" ./scatter_sfr_block_linear_$i.in
  sed -i "19s/scatter_sfr_block_linear_0/scatter_sfr_block_linear_$i/g" ./scatter_sfr_block_linear_$i.in
done