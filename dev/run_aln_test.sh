#!/usr/bin/env bash
echo "Running kalign basic alignment tests:";


testfiles=("BB11001.msf" "BB11001_EOF.msf" "BB11001.tfa" "BB12006.msf" "BB12006.tfa" "BB30014.msf" "BB30014.tfa")

for i in "${testfiles[@]}"
do
    echo $i

    error=$( ../src/kalign --devtest $i -o test.afa 2>&1 )
    status=$?
    if [[ $status -eq 0 ]]; then
        printf "%10s%10s%10s\n"  read/write $i SUCCESS;
    else
        printf "%10s%10s%10s\n"  read/write $i FAILED;
        printf "with ERROR $status and Message:\n\n$error\n\n";
        exit 1;
    fi

done
