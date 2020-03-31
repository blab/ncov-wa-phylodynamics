#!/bin/bash

for clock in 8e-4; do

    for bu in 36.5 30.4; do

        # for idx in `seq 1 5`; do
        for idx in `seq 1 3`; do

            java -jar ~/code/beast_and_friends/EpiInf/out/artifacts/EpiInf_jar/EpiInf.jar \
                 -overwrite \
                 -D outbreak=$outbreak,clock=$clock,bu=$bu,idx=$idx \
                 TrajectoryMapper.xml

        done
        
    done
    
done
