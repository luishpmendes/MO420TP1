#!/bin/bash
path="../..";
make outputAggregator;
echo "Aggregating outputs";
for k in 1 2
do
    echo "k = "$k;
    ./outputAggregator $path"/outputs/"$k"/type1/z50-200-199.out" $path"/outputs/"$k"/type1/z50-200-398.out" $path"/outputs/"$k"/type1/z50-200-597.out" $path"/outputs/"$k"/type1/z50-200-995.out" $path"/outputs/"$k"/type1/z100-300-448.out" $path"/outputs/"$k"/type1/z100-300-897.out" $path"/outputs/"$k"/type1/z100-500-1247.out" $path"/outputs/"$k"/type1/z100-500-2495.out" $path"/outputs/"$k"/type1/z100-500-3741.out" $path"/outputs/"$k"/type1/z200-600-1797.out" $path"/outputs/"$k"/type1/z200-800-3196.out" $path"/outputs/"$k"/type2/z50-200-3903.out" $path"/outputs/"$k"/type2/z50-200-4877.out" $path"/outputs/"$k"/type2/z50-200-5864.out" $path"/outputs/"$k"/type2/z100-300-8609.out" $path"/outputs/"$k"/type2/z100-300-10686.out" $path"/outputs/"$k"/type2/z100-300-12761.out" $path"/outputs/"$k"/type2/z100-500-24740.out" $path"/outputs/"$k"/type2/z100-500-30886.out" $path"/outputs/"$k"/type2/z100-500-36827.out" $path"/outputs/"$k"/type2/z200-400-13660.out" $path"/outputs/"$k"/type2/z200-400-17089.out" $path"/outputs/"$k"/type2/z200-400-20469.out" $path"/outputs/"$k"/type2/z200-600-34504.out" $path"/outputs/"$k"/type2/z200-600-42860.out" $path"/outputs/"$k"/type2/z200-600-50984.out" $path"/outputs/"$k"/type2/z200-800-62625.out" $path"/outputs/"$k"/type2/z200-800-78387.out" $path"/outputs/"$k"/type2/z200-800-93978.out" $path"/outputs/"$k"/type2/z300-600-31000.out" $path"/outputs/"$k"/type2/z300-600-38216.out" $path"/outputs/"$k"/type2/z300-600-45310.out" $path"/outputs/"$k"/type2/z300-800-59600.out" $path"/outputs/"$k"/type2/z300-800-74500.out" $path"/outputs/"$k"/type2/z300-800-89300.out" $path"/outputs/"$k"/type2/z300-1000-96590.out" $path"/outputs/"$k"/type2/z300-1000-120500.out" $path"/outputs/"$k"/type2/z300-1000-144090.out" > $path"/outputs/"$k"/relaxlag.out";
done
make clean;

