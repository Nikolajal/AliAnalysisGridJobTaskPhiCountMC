#! /bin/bash

echo "starting: $@" > starter.log
echo "start time: $(date)" >> starter.log
$@
echo "end time: $(date)" >> starter.log
echo "workdir:" >> starter.log
ls -l >> starter.log
