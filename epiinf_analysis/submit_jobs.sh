#!/bin/bash

job_prefix=report3003

if ! [ -d Results ]; then
    mkdir Results
fi
for clock in 8e-4; do

    for bu in 36.5 30.4; do

        jobid="${job_prefix}_${clock}_${bu}_BDSKY"

        bsub <<EOF
#!/bin/sh
#BSUB -W 24:00
#BSUB -R "rusage[mem=4096]"
#BSUB -J "$jobid[1-5]"
module load java
JAVA="java -Xmx3G"
JAR=\$HOME/bdsky.jar

SEED=\$LSB_JOBINDEX
STATEFILE=Results/BD_epi.clock_$clock.bu_$bu.\$SEED.state

\$JAVA -jar \$JAR -seed \$SEED -D clockrate=$clock,burate=$bu -statefile \$STATEFILE -overwrite BD_epi.xml

EOF

    done

done


