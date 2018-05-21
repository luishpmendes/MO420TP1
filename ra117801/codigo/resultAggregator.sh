#!/bin/bash
path="../../";
make resultAggregator;
echo "Aggregating results";
for k in 1 2
do
    echo "k = "$k;
    ./resultAggregator $path"/results/"$k"/type1/z50-200-199.res" $path"/results/"$k"/type1/z50-200-398.res" $path"/results/"$k"/type1/z50-200-597.res" $path"/results/"$k"/type1/z50-200-995.res" $path"/results/"$k"/type1/z100-300-448.res" $path"/results/"$k"/type1/z100-300-897.res" $path"/results/"$k"/type1/z100-500-1247.res" $path"/results/"$k"/type1/z100-500-2495.res" $path"/results/"$k"/type1/z100-500-3741.res" $path"/results/"$k"/type1/z200-600-1797.res" $path"/results/"$k"/type1/z200-800-3196.res" $path"/results/"$k"/type2/z50-200-3903.res" $path"/results/"$k"/type2/z50-200-4877.res" $path"/results/"$k"/type2/z50-200-5864.res" $path"/results/"$k"/type2/z100-300-8609.res" $path"/results/"$k"/type2/z100-300-10686.res" $path"/results/"$k"/type2/z100-300-12761.res" $path"/results/"$k"/type2/z100-500-24740.res" $path"/results/"$k"/type2/z100-500-30886.res" $path"/results/"$k"/type2/z100-500-36827.res" $path"/results/"$k"/type2/z200-400-13660.res" $path"/results/"$k"/type2/z200-400-17089.res" $path"/results/"$k"/type2/z200-400-20469.res" $path"/results/"$k"/type2/z200-600-34504.res" $path"/results/"$k"/type2/z200-600-42860.res" $path"/results/"$k"/type2/z200-600-50984.res" $path"/results/"$k"/type2/z200-800-62625.res" $path"/results/"$k"/type2/z200-800-78387.res" $path"/results/"$k"/type2/z200-800-93978.res" $path"/results/"$k"/type2/z300-600-31000.res" $path"/results/"$k"/type2/z300-600-38216.res" $path"/results/"$k"/type2/z300-600-45310.res" $path"/results/"$k"/type2/z300-800-59600.res" $path"/results/"$k"/type2/z300-800-74500.res" $path"/results/"$k"/type2/z300-800-89300.res" $path"/results/"$k"/type2/z300-1000-96590.res" $path"/results/"$k"/type2/z300-1000-120500.res" $path"/results/"$k"/type2/z300-1000-144090.res" > $path"/results/"$k"/relaxlag.res";
done
make clean;

