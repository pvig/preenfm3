# preenfm3

Fork from https://github.com/Ixox/preenfm3  

## Some Reverb ~ Strange 8

##### the Dattoro Reverb blocks  
  
input  
quad allpass filter  
lp & hp  
pre delay  
input diffuser   
reverb tank ( (allpass>delay>lp>allpass>delay) > (allpass>delay>allpass>delay) )  
output from various point in the tank  

--------------

##### Mix page 3 : Dry/wet

Set fx level for each instrument  

--------------

##### Global page 3 : Master Fx 1

* Preset
    * Group of parameters :
        * Size
        * Lfo speed
        * Lfo Depth
        * Decay
        * Diffusion
        * Notch Base
        * Damping
* Pre delay 
    * Time of pre delay
* Pre delay Mix 
    * the input diffuser is fed with a mix of audio input and pre delayed input.
* Filter Base
    * Input filter hp/lp
* Filter Width 
* Notch Spread 
    * Input is filtered with a 4 pole allpass filter

--------------

##### Midi CC : on each timbre midi channel

11.    Send level

##### Midi CC : on Global midi channel

40.    PRESET
41.    PREDELAYTIME
42.    PREDELAY MIX
43.    TILT
44.    LFO SPEED
45.    LFO DEPTH

--------------

## Credits
JON DATTORRO for the reverb algorithm
https://ccrma.stanford.edu/~dattorro/EffectDesignPart1.pdf

Dale Johnson of ValleyAudio for some great code example & ideas from the Plateau reverb
https://github.com/ValleyAudio/ValleyRackFree/  
