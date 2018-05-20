#!/bin/bash
make relaxlag;
path="../../";
for k in 1 2
do
    echo "k = "$k;
    mkdir -p $path"results/"$k"/";
    mkdir -p $path"results/"$k"/type1/";
    mkdir -p $path"results/"$k"/type2/";
    mkdir -p $path"outputs/"$k"/";
    mkdir -p $path"outputs/"$k"/type1/";
    mkdir -p $path"outputs/"$k"/type2/";
    for instance in "type1/z50-200-199" "type1/z50-200-398" "type1/z50-200-597" "type1/z50-200-995" "type1/z100-300-448" "type1/z100-300-897" "type1/z100-500-1247" "type1/z100-500-2495" "type1/z100-500-3741" "type1/z200-600-1797" "type1/z200-800-3196" "type2/z50-200-3903" "type2/z50-200-4877" "type2/z50-200-5864" "type2/z100-300-8609" "type2/z100-300-10686" "type2/z100-300-12761" "type2/z100-500-24740" "type2/z100-500-30886" "type2/z100-500-36827" "type2/z200-400-13660" "type2/z200-400-17089" "type2/z200-400-20469" "type2/z200-600-34504" "type2/z200-600-42860" "type2/z200-600-50984" "type2/z200-800-62625" "type2/z200-800-78387" "type2/z200-800-93978" "type2/z300-600-31000" "type2/z300-600-38216" "type2/z300-600-45310" "type2/z300-800-59600" "type2/z300-800-74500" "type2/z300-800-89300" "type2/z300-1000-96590" "type2/z300-1000-120500" "type2/z300-1000-144090"
    do
        echo "instance = "$instance;
        ./relaxlag $k $path"instances/"$instance".gcc" $path"/results/"$k"/"$instance".res" > $path"/outputs/"$k"/"$instance".out";
    done
    echo "Aggregating results";
    make resultAggregator;
    ./resultAggregator $path"/results/"$k"/type1/z50-200-199.res" $path"/results/"$k"/type1/z50-200-398.res" $path"/results/"$k"/type1/z50-200-597.res" $path"/results/"$k"/type1/z50-200-995.res" $path"/results/"$k"/type1/z100-300-448.res" $path"/results/"$k"/type1/z100-300-897.res" $path"/results/"$k"/type1/z100-500-1247.res" $path"/results/"$k"/type1/z100-500-2495.res" $path"/results/"$k"/type1/z100-500-3741.res" $path"/results/"$k"/type1/z200-600-1797.res" $path"/results/"$k"/type1/z200-800-3196.res" $path"/results/"$k"/type2/z50-200-3903.res" $path"/results/"$k"/type2/z50-200-4877.res" $path"/results/"$k"/type2/z50-200-5864.res" $path"/results/"$k"/type2/z100-300-8609.res" $path"/results/"$k"/type2/z100-300-10686.res" $path"/results/"$k"/type2/z100-300-12761.res" $path"/results/"$k"/type2/z100-500-24740.res" $path"/results/"$k"/type2/z100-500-30886.res" $path"/results/"$k"/type2/z100-500-36827.res" $path"/results/"$k"/type2/z200-400-13660.res" $path"/results/"$k"/type2/z200-400-17089.res" $path"/results/"$k"/type2/z200-400-20469.res" $path"/results/"$k"/type2/z200-600-34504.res" $path"/results/"$k"/type2/z200-600-42860.res" $path"/results/"$k"/type2/z200-600-50984.res" $path"/results/"$k"/type2/z200-800-62625.res" $path"/results/"$k"/type2/z200-800-78387.res" $path"/results/"$k"/type2/z200-800-93978.res" $path"/results/"$k"/type2/z300-600-31000.res" $path"/results/"$k"/type2/z300-600-38216.res" $path"/results/"$k"/type2/z300-600-45310.res" $path"/results/"$k"/type2/z300-800-59600.res" $path"/results/"$k"/type2/z300-800-74500.res" $path"/results/"$k"/type2/z300-800-89300.res" $path"/results/"$k"/type2/z300-1000-96590.res" $path"/results/"$k"/type2/z300-1000-120500.res" $path"/results/"$k"/type2/z300-1000-144090.res" > $path"/results/"$k"/relaxlag.res";
done
make clean;

