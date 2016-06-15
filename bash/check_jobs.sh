#!/bin/bash

# TODO-LIST: READ *ONLY* THE SIM NUMBER FROM GENERAL SCRIPT NAMES

IC2CL=xt8786@ic2.scc.kit.edu
COMMON_PATH=/work/kit/istm/xt8786/EFLUX/
EMAIL=davide.gatti@kit.edu
JOB_NAME[0]="CPI-OW-1-2";   JOB_SCRIPT[0]="CPI/RePI-6500/OW/script-1-2.sh";
JOB_NAME[1]="CPI-OW-3-4";   JOB_SCRIPT[1]="CPI/RePI-6500/OW/script-3-4.sh";
JOB_NAME[2]="CPI-OW-5-6";   JOB_SCRIPT[2]="CPI/RePI-6500/OW/script-5-6.sh";
JOB_NAME[3]="CPI-OW-7-8";   JOB_SCRIPT[3]="CPI/RePI-6500/OW/script-7-8.sh";
JOB_NAME[4]="CPI-OW-9-10";  JOB_SCRIPT[4]="CPI/RePI-6500/OW/script-9-10.sh";
JOB_NAME[5]="CPI-VC-1-2";   JOB_SCRIPT[5]="CPI/RePI-6500/VC/script-1-2.sh";
JOB_NAME[6]="CPI-VC-3-4";   JOB_SCRIPT[6]="CPI/RePI-6500/VC/script-3-4.sh";
JOB_NAME[7]="CPI-VC-5-6";   JOB_SCRIPT[7]="CPI/RePI-6500/VC/script-5-6.sh";
JOB_NAME[8]="CtPI-VC-1-2";  JOB_SCRIPT[8]="CtPI/RePI-6500/VC/script-1-2.sh";
JOB_NAME[9]="CtPI-VC-3-4";  JOB_SCRIPT[9]="CtPI/RePI-6500/VC/script-3-4.sh";
JOB_NAME[10]="CtPI-VC-5-6"; JOB_SCRIPT[10]="CtPI/RePI-6500/VC/script-5-6.sh";
JOB_NAME[11]="CPI-TW-1-2";  JOB_SCRIPT[11]="CPI/RePI-6500/TW/script-1-2.sh";
JOB_NAME[12]="CPI-TW-3-4";  JOB_SCRIPT[12]="CPI/RePI-6500/TW/script-3-4.sh";
JOB_NAME[13]="CPI-TW-5-6";  JOB_SCRIPT[13]="CPI/RePI-6500/TW/script-5-6.sh";

freq=1200

RED=$(tput bold; tput setaf 1)
GREEN=$(tput bold; tput setaf 2)
CYAN=$(tput bold; tput setaf 6)
NORM=$(tput sgr0)

while true; do

  clear; echo ""
  tput bold; tput setaf 2; echo "=====> Cheching the status of the jobs on IC2"; tput sgr0
  printf "%15s %20s\n" "Job name" "Status"
  printf "%15s %20s\n" "--------" "------"
  out=$( ssh -tq $IC2CL '/jms/bin/job_queue -l' )
  for cases in {0..13}; do
        job=$( echo "$out" | grep -w "${JOB_NAME[$cases]}" | awk '{print $4}' )
        stat=$( echo "$out" | grep -w "${JOB_NAME[$cases]}" |  awk '{print $14}' )
        if [ "$job" == "" ]; then
          isfinished=yes
          casel=`echo ${JOB_SCRIPT[$cases]} | sed 's/[^0-9]/ /g' | awk '{print $2}'`;
          caseh=`echo ${JOB_SCRIPT[$cases]} | sed 's/[^0-9]/ /g' | awk '{print $3}'`;
          for ((i=casel; i<=caseh; i++)); do
              simpath="$COMMON_PATH$(echo ${JOB_SCRIPT[$cases]} | sed -e 's/\(.*\)script.*/\1/')$i"
              t_sim=$(ssh $IC2CL "tail -1 ${simpath}/Runtimedata" | awk '{print $1}');
              t_max=$(ssh $IC2CL "grep -w t_max ${simpath}/dns.in" | sed -e 's/.*t_max=\(.*\)time.*/\1/');
              dt=$(ssh $IC2CL "tail -1 ${simpath}/Runtimedata" | awk '{print $11}');
              if [ "$t_sim" == ""  ]; then
                 isfinished=no;
              else
                 if [ $(echo "$t_sim"'>'"$t_max"'-'"$dt"'/2' | bc -l) == 0 ]; then
                 isfinished=no;
                 fi
              fi
          done
          if [ "$isfinished" == "yes" ]; then
              stat="${GREEN}Finished!${NORM}"
              #ssh $IC2CL "mailx -s \"Case ${JOB_NAME[$cases]} FINISHED on IC2\" < /dev/null $EMAIL"
          else
              stat="${RED}Not in queue:${NORM} restarting"
              ssh $IC2CL "mailx -s \"Case ${JOB_NAME[$cases]} is not running on IC2: restarting\" < /dev/null $EMAIL"
              dir="$COMMON_PATH$(echo ${JOB_SCRIPT[$cases]} | sed -e 's/\(.*\)script.*/\1/')";
	      ssh $IC2CL "cd ${dir}; ${COMMON_PATH}${JOB_SCRIPT[$cases]}" > /dev/null &
          fi

        fi
        printf "%15s %20s\n" "${JOB_NAME[$cases]}" "$stat"
  done


echo ""
echo "${GREEN}====> Job scan completed, I will retry automatically in ${CYAN}$((freq/60))${GREEN} minutes... ${NORM}"
sleep $freq
done
