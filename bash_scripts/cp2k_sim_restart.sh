#!/bin/bash

#cp PYRENE_1H2O-pos-1.xyz pos.xyz; cp PYRENE_1H2O-vel-1.xyz vel.xyz

VEL_POS_LEN=31
DIP_LEN=6
CHARGE_LEN=37


LAST_STEP=$(tac PYRENE_1H2O-pos-1.xyz | grep -m 1 'i =' | awk -F '[^0-9]*' '{print $2}')

FIRST_STEP=$(grep -m 1 'i =' PYRENE_1H2O-pos-1.xyz | awk -F '[^0-9]*' '{print $2}')

START_STEP=$(awk '/STEP_START_VAL/ {print $2}' PYRENE_1H2O-1.restart)


if [[ $LAST_STEP == $START_STEP ]]; then
	tail -n +$(($DIP_LEN + 1)) dipole.txt > dip.txt
	tail -n +$((CHARGE_LEN + 1)) hirshfeld.txt > hs.txt
	cp PYRENE_1H2O-pos-1.xyz pos.xyz
	cp PYRENE_1H2O-vel-1.xyz vel.xyz
	#cp dipole.txt dip.txt
	#cp hirshfeld.txt hs.txt

else
	DEL_DIFF=$(($LAST_STEP - $START_STEP))
	DIP_REM=$(($DIP_LEN*$DEL_DIFF))
	VEL_POS_REM=$(($VEL_POS_LEN*$DEL_DIFF))
	CHARGE_REM=$(($CHARGE_LEN*$DEL_DIFF))
	tail -n +$(($DIP_LEN + 1)) dipole.txt |  head -n -$DIP_REM > dip.txt
	tail -n +$((CHARGE_LEN + 1)) hirshfeld.txt |head -n -$CHARGE_REM  > hs.txt
	head -n -$VEL_POS_REM PYRENE_1H2O-pos-1.xyz > pos.xyz
	head -n -$VEL_POS_REM PYRENE_1H2O-vel-1.xyz > vel.xyz
	#head -n -$DIP_REM dipole.txt > dip.txt
	#head -n -$CHARGE_REM hirshfeld.txt  > hs.txt
fi

echo $LAST_STEP
echo $(($START_STEP - $FIRST_STEP + 1))

HS_LINES=$(wc -l < hs.txt)
HS_DIFF=$((HS_LINES / 37))

if [[ $(($START_STEP - $FIRST_STEP + 1)) == $HS_DIFF ]]; then
	echo "All Good!"

else
	echo "Oh no! There's a problem!"
	echo $(($START_STEP - $FIRST_STEP + 1))
	echo $HS_DIFF
fi

#mkdir next

#cp {pyr_start.txt,md_pyrene_1w.inp,PYRENE_1H2O-RESTART.wfn,PYRENE_1H2O-1.restart,run_cp2k.sh} next/


