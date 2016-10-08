#~/bin/sh

for file in $(ls code)
do
    #this will work on the VLSCI clusters, if on OSX you see the error

    #./code/base.cpp:217:21: error: expected expression
    #    cities[i] = {x, y};

    # then brew install gcc-6
    # and replace g++ with /usr/local/bin/g++-6

    /usr/local/bin/g++-6 ./code/$file -o ./execs/${file%.*}
done
