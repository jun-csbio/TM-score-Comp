#!/bin/bash

# If the following compile is failed, please try:
#g++ -o ./exe/TMscoreC ./code/TMscoreC.cpp
g++ -static -O3 -ffast-math -lm -o ./exe/TMscoreC ./code/TMscoreC.cpp

# If you want to recompile the multi-thread TM-score-Comp, please run
# the following command. 
#g++ -static -O3 -ffast-math -lm -o ./exe/TMscoreCmt ./code/TMscoreCmt.cpp
# Or run the command of
#g++ -o ./exe/TMscoreCmt ./code/TMscoreCmt.cpp
# If you can not compile them, you should update your g++ version or 
# directly use the executable of ./exe/TMscoreCmt

chmod -R 777 ./exe/
