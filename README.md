# preenfm3

Fork from https://github.com/Ixox/preenfm3  

## Some Reverb ~ Strange 8

##### the Dattoro Reverb blocks  
  
input  
lp & hp (tilt)  
pre delay  
input diffuser   
reverb tank ( (allpass>delay>lp>allpass>delay) > (allpass>delay>allpass>delay) )  
output from various point in the tank  

--------------

##### Mix page 3 : Dry/wet

Set fx level for each instrument  

--------------

##### Global page 3 : Master Fx 1

* Pre delay 
    * Time of pre delay
* Pre delay Mix 
    * the input diffuser is fed with a mix of audio input and pre delayed input.
* Diffusion 
    * set the input & tank diffuser allpass coefficient
* Decay 
    * tank feedback
* Size 
    * size of the tank delay
* Damping 
    * cut off of the tank lowpass

--------------

##### Global page 4 : Master Fx 2

* Mod Speed 
    * lfo speed (4 of them)
* Mod Depth 
    * lfo mod depth, modulate the 4 delay allpass in the tank
* Env Mod 
    * envelope mod depth, make the tail pitch shifting up or down
* Threshold 
    * audio level at which envelope trigger
* Input Tilt 
    * input filter from lowpass ~0 to high pass ~1
* Env Feedback 
    * feedback (decay) modulation with the envelope

the envelope release rate is relative to the decay param.  

--------------

##### Midi CC : on Global midi channel

34.    Send level 1
35.    Send level 2
36.    Send level 3
37.    Send level 4
38.    Send level 5
39.    Send level 6

40.    not sure yet..
41.    ,,,
42.    ,,,
43.    ,,,
44.    ,,,
45.    ,,,   

--------------

## Credits
JON DATTORRO for the reverb algorithm
https://ccrma.stanford.edu/~dattorro/EffectDesignPart1.pdf

Dale Johnson of ValleyAudio for lot of code inspiration from the Plateau reverb
https://github.com/ValleyAudio/ValleyRackFree/