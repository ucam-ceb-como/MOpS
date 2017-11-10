#!/bin/bash

nsp=$1
m0=$2
nruns=$3
nparts=$4

inp="mops-hm-s1.xml"
tmp1="mops-hm-temp-dz.xml"
tmp2="mops-hm-temp-wz.xml"
tmp3="mops-hm-temp-cz.xml"

tmp0="tmops-hm-temp.xml"


# input file
cp $inp $tmp0

eval "sed 4's/.*/    <runs>$nruns<\/runs>/' $inp" > $tmp0	
eval "sed ' ' $tmp0" > "$inp"
eval "sed 8's/.*/    <pcount>$nsp<\/pcount>/' $inp" > $tmp0	
eval "sed ' ' $tmp0" > "$inp"
eval "sed 9's/.*/    <maxm0>$m0<\/maxm0>/' $inp" > $tmp0	
eval "sed ' ' $tmp0" > "$inp"
eval "sed 63's/.*/          <ptrack enable=\"true\" ptcount=\"$nparts\"\/>/' $inp" > $tmp0	
eval "sed ' ' $tmp0" > "$inp"

rm -f $tmp0


# DZ temp file
cp $tmp1 $tmp0

eval "sed 4's/.*/    <runs>$nruns<\/runs>/' $tmp1" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp1"
eval "sed 8's/.*/    <pcount>$nsp<\/pcount>/' $tmp1" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp1"
eval "sed 9's/.*/    <maxm0>$m0<\/maxm0>/' $tmp1" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp1"
eval "sed 112's/.*/        <ptrack enable=\"true\" ptcount=\"$nparts\"\/>/' $tmp1" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp1"

rm -f $tmp0


# WZ temp file
cp $tmp2 $tmp0

eval "sed 3's/.*/    <runs>$nruns<\/runs>/' $tmp2" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp2"
eval "sed 7's/.*/    <pcount>$nsp<\/pcount>/' $tmp2" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp2"
eval "sed 8's/.*/    <maxm0>$m0<\/maxm0>/' $tmp2" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp2"
eval "sed 69's/.*/        <ptrack enable=\"true\" ptcount=\"$nparts\"\/>/' $tmp2" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp2"

rm -f $tmp0


# CZ temp file
cp $tmp3 $tmp0

eval "sed 3's/.*/    <runs>$nruns<\/runs>/' $tmp3" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp3"
eval "sed 7's/.*/    <pcount>$nsp<\/pcount>/' $tmp3" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp3"
eval "sed 8's/.*/    <maxm0>$m0<\/maxm0>/' $tmp3" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp3"
eval "sed 73's/.*/        <ptrack enable=\"true\" ptcount=\"$nparts\"\/>/' $tmp3" > $tmp0	
eval "sed ' ' $tmp0" > "$tmp3"

rm -f $tmp0








