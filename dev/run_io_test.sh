#!/usr/bin/env bash
echo "Running kalign basic io test:";


testfiles=("a2m.good.1" "a2m.good.2" "afa.good.1" "afa.good.2" "afa.good.3" "clustal.good.1" "clustal.good.2")


for i in "${testfiles[@]}"
do
    echo $i

    error=$( ../src/rwaln $i   2>&1 )
    status=$?
    if [[ $status -eq 0 ]]; then
        printf "%10s%10s%10s\n"  read/write $i SUCCESS;
    else
        printf "%10s%10s%10s\n"  read/write $i FAILED;
        printf "with ERROR $status and Message:\n\n$error\n\n";
        exit 1;
    fi
done
