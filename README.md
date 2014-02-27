HarmOscillChain
===============
In the block analysis the program receive as input the frequencies and seeds files.

Update: In the all in one versione the main receive this argument in THIS order:
1)frequencies vector file
2)total simulation time
3)type of themrostat to use
4)thermostat relaxation time (tau)
5)a boolean variable: 1=create a file with all the trajectory formatted as: q1 p1 q2
6) number of runs to compute if block analysis is requested
7)file with seeds to use in block analysis

5) 6) 7) are not necessary if block analysis is not needed!
