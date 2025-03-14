#!/bin/bash

echo First run the compilation with
echo make clean
echo 'make -j`nproc` 2>&1 | tee build.log'
echo ''
echo 'Current results:'
egrep '(warning:|error:)' build.log | sed -e 's/: .*\[/ [/' | sort -u | sed -e 's/.*\[/[/g' | sort | uniq -c
echo ''
echo 'In detail:'
egrep '(warning:|error:)' build.log | sed -e 's/: .*\[/ [/' | sort -u
