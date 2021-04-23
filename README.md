# preenfm3

Fork from https://github.com/Ixox/preenfm3


## master fx branch


#### Inspired by the pure data/cyclone comb filter :



![Alt text](/doc/comb_cyclone.png?raw=true "Comb filter from pure data/cyclone")


(screenshot from the help file of the cyclone object)


#### input of the comb filter is fed with an harmonic tremolo



### The filter is controlable with 3 menu pages


#### Mix page : Fx send

Set send level for each instrument


#### Global page : Master Fx 1

Comb filter parameters (see screenshot)

1. Time
2. FeedForward
3. Feedback
4. Input
5. Threshold
6. Release

Threshold & Release params are for the envelope follower
If audio is over threshold, the envelope reach its maximum very fast.
Then, when signal go below threshold, the envelope enter release state.

#### Global page : Master Fx 2

Modulation parameters

1. FeedBack Speed 
2. FeedBack Depth 
3. FeedBack Env mod
4. Tremolo Speed 
5. Tremolo Depth 
6. Tremolo Env mod


##### Midi CC control change

34.    SEND1
35.    SEND2
36.    SEND3
37.    SEND4
38.    SEND5
39.    SEND6
40.    TIME
41.    FEEDFORWARD
42.    FEEDBACK
43.    INPUT LEVEL
44.    MOD LEVEL
45.    MOD SPEED   

