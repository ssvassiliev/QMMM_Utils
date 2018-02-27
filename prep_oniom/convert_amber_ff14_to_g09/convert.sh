#! /bin/bash
# Convert parm10.dat
# VDW
echo "!" > ff14_g09.prm
echo "NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.000" >> ff14_g09.prm
echo "!" >> ff14_g09.prm
echo "!" >> ff14_g09.prm
echo "! Vanderwaals parameters" >> ff14_g09.prm
echo "!" >> ff14_g09.prm
awk 'NR==962,NR==1000' parm10.dat | awk '{printf "VDW  " "%-4s %-10s %-10s\n",$1,$2,$3}' >> ff14_g09.prm
# BONDS
echo "!" >> ff14_g09.prm
echo "! Stretches" >> ff14_g09.prm
echo "!" >> ff14_g09.prm
 awk 'NR>=67 && NR<=217 {gsub(/-/," "); printf "HrmStr1  " "%-4s %-4s %-10s %-10s\n",$1,$2,$3,$4}' parm10.dat >> ff14_g09.prm
# ANGLES
echo "!" >> ff14_g09.prm
echo "! Bendings" >> ff14_g09.prm
echo "!" >> ff14_g09.prm
 awk 'NR>=219 && NR<=618 {gsub(/-/," "); printf "HrmBnd1  " "%-4s %-4s %-4s %-8s %-10s\n",$1,$2,$3,$4,$5}' parm10.dat >> ff14_g09.prm
# TORSIONS 
awk 'NR>=620 && NR<=894' parm10.dat > torsions.dat
gcc ambtrs2gau.c -o AmbTrs2Gau
./AmbTrs2Gau >> ff14_g09.prm
rm torsions.dat AmbTrs2Gau
# IMPROPER TORSIONS
echo "!" >> ff14_g09.prm
echo "! Improper torsions" >> ff14_g09.prm
echo "!" >> ff14_g09.prm
awk 'NR>=896 && NR<=954 {gsub(/-/," "); gsub(/X/,"*"); printf "ImpTrs  " "%-4s %-4s %-4s %-4s %6.1f %8.1f %8.1f\n",$1,$2,$3,$4,$5,$6,$7}' parm10.dat >> ff14_g09.prm
#
# Convert frcmod.ff14SB
# VDW
echo "!" > ff14SB_g09.prm
echo "! Vanderwaals parameters" >> ff14SB_g09.prm
echo "!" >> ff14SB_g09.prm
awk 'NR==509,NR==512' frcmod.ff14SB | awk '{printf "VDW  " "%-6s %-10s %-10s\n",$1,$2,$3}' >> ff14SB_g09.prm
# BONDS
echo "!" >> ff14SB_g09.prm
echo "! Stretches" >> ff14SB_g09.prm
echo "!" >> ff14SB_g09.prm
 awk 'NR>=9 && NR<=35 {gsub(/-/," "); printf "HrmStr1  " "%-6s %-4s %-10s %-10s\n",$1,$2,$3,$4}' frcmod.ff14SB >> ff14SB_g09.prm
# ANGLE
echo "!" >> ff14SB_g09.prm
echo "! Bendings" >> ff14SB_g09.prm
echo "!" >> ff14SB_g09.prm
 awk 'NR>=38 && NR<=130 {gsub(/-/," "); printf "HrmBnd1  " "%-6s %-4s %-4s %-8s %-10s\n",$1,$2,$3,$4,$5}' frcmod.ff14SB >> ff14SB_g09.prm
# TORSIONS 
awk 'NR>=133 && NR<=501' frcmod.ff14SB > torsions.dat
gcc ambtrs2gau.c -o AmbTrs2Gau
./AmbTrs2Gau >> ff14SB_g09.prm
rm torsions.dat AmbTrs2Gau
# IMPROPER TORSIONS
echo "!" >> frcmod-can.prm
echo "! Improper torsions" >> frcmod-can.prm
echo "!" >> frcmod-can.prm
awk 'NR>=504 && NR<=506 {gsub(/-/," "); gsub(/X/,"*"); printf "ImpTrs  " "%-6s %-4s %-4s %-4s %6.1f %8.1f %8.1f\n",$1,$2,$3,$4,$5,$6,$7}' frcmod-can.dat >> frcmod-can.prm
#
# Convert frcmod-can-01.dat
# VDW
echo "!" > frcmod-can.prm
echo "! Vanderwaals parameters" >> frcmod-can.prm
echo "!" >> frcmod-can.prm
awk 'NR==19,NR==31' frcmod-can.dat | awk '{printf "VDW  " "%-6s %-10s %-10s\n",$1,$2,$3}' >> frcmod-can.prm
# BONDS
echo "!" >> frcmod-can.prm
echo "! Stretches" >> frcmod-can.prm
echo "!" >> frcmod-can.prm
 awk 'NR>=34 && NR<=63 {gsub(/-/," "); printf "HrmStr1  " "%-6s %-4s %-10s %-10s\n",$1,$2,$3,$4}' frcmod-can.dat >> frcmod-can.prm
# ANGLE
echo "!" >> frcmod-can.prm
echo "! Bendings" >> frcmod-can.prm
echo "!" >> frcmod-can.prm
 awk 'NR>=66 && NR<=130 {gsub(/-/," "); printf "HrmBnd1  " "%-6s %-4s %-4s %-8s %-10s\n",$1,$2,$3,$4,$5}' frcmod-can.dat >> frcmod-can.prm
# TORSIONS 
awk 'NR>=133 && NR<=163' frcmod-can.dat > torsions.dat
gcc ambtrs2gau.c -o AmbTrs2Gau
./AmbTrs2Gau >> frcmod-can.prm
rm torsions.dat AmbTrs2Gau
# IMPROPER TORSIONS
echo "!" >> frcmod-can.prm
echo "! Improper torsions" >> frcmod-can.prm
echo "!" >> frcmod-can.prm
awk 'NR>=166 && NR<=167 {gsub(/-/," "); gsub(/X/,"*"); printf "ImpTrs  " "%-6s %-4s %-4s %-4s %6.1f %8.1f %8.1f\n",$1,$2,$3,$4,$6,$7,$5}' frcmod-can.dat >> frcmod-can.prm











