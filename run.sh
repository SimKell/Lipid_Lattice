#!/bin/bash

FILE="input.txt"

Conc_0="(CHOL): 0.0"
Conc_1="(CHOL): 0.1"
Conc_2="(CHOL): 0.2"
Conc_3="(CHOL): 0.3"

Temp_0="Temperature: 310"
Temp_1="Temperature: 300"
Temp_2="Temperature: 290"

Proj_0="Name: Test"
Proj_1="Name: 00_310"
Proj_2="Name: 00_300"
Proj_3="Name: 00_290"

sed -i "s/$Proj_0/$Proj_1/" "$FILE"

./a.out
sed -i "s/$Conc_0/$Conc_1/" "$FILE"
mv "00_310/0" "00_310/00"

./a.out
sed -i "s/$Conc_1/$Conc_2/" "$FILE"
mv "00_310/0" "00_310/01"

./a.out
sed -i "s/$Conc_2/$Conc_3/" "$FILE"
mv "00_310/0" "00_310/02"

./a.out
sed -i "s/$Conc_3/$Conc_0/" "$FILE"
mv "00_310/0" "00_310/03"



sed -i "s/$Proj_1/$Proj_2/" "$FILE"
sed -i "s/$Temp_0/$Temp_1/" "$FILE"

./a.out
sed -i "s/$Conc_0/$Conc_1/" "$FILE"
mv "00_300/0" "00_300/00"

./a.out
sed -i "s/$Conc_1/$Conc_2/" "$FILE"
mv "00_300/0" "00_300/01"

./a.out
sed -i "s/$Conc_2/$Conc_3/" "$FILE"
mv "00_300/0" "00_300/02"

./a.out
sed -i "s/$Conc_3/$Conc_0/" "$FILE"
mv "00_300/0" "00_300/03"



sed -i "s/$Proj_2/$Proj_3/" "$FILE"
sed -i "s/$Temp_1/$Temp_2/" "$FILE"

./a.out
sed -i "s/$Conc_0/$Conc_1/" "$FILE"
mv "00_290/0" "00_290/00"

./a.out
sed -i "s/$Conc_1/$Conc_2/" "$FILE"
mv "00_290/0" "00_290/01"

./a.out
sed -i "s/$Conc_2/$Conc_3/" "$FILE"
mv "00_290/0" "00_290/02"

./a.out
sed -i "s/$Conc_3/$Conc_0/" "$FILE"
mv "00_290/0" "00_290/03"