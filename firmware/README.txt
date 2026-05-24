TODO

- Clean old editor completely.
- Put back audio engine reinit (BACK + MENU).
- Move Gate FX to regular filters.

Bug

- None listed yet.

In Progress

- Finish all menus.

Not Tested

- None listed yet.

-----------------------------------------------------
Optional / Nice To Have

- UI: show LFO/OSC shapes.
- UI: show env (other envs?).
- MIDI send and receive per instrument in mixer.
- Editor OP: display OP number in bottom middle empty button (e.g. 1/3).
-----------------------------------------------------

EPIC

- Add feedback to all algos?
- Add limiter for each stereo DAC.

Done

- 2 destinations per matrix row.
- CPU usage as settings option + TFT access outside audio thread.
- Oscilloscope moves when fTines != 0. Use real frequency in Osc.cpp.
- Add all/current instruments MIDI input in mixer/other menu.
- Get rid of MAX_NUMBER_OF_OPERATORS. Use MAX_NUMBER_OF_VOICES only.
- Update screen when new value comes from external input.
- Bug: Reset ('-' + encoder) does not work on OP envelope (for OP2 and above).
- Bug: changing level on OP2+ calls reloadASDR for OP1.env only.
- Mono(2) plays 2 voices. It should play only 1.
- Display OP oscillator Ftyp = fixed.
- Display LFO frequency: after 99.5 / M/16 etc., erase problem fixed.
- Display LFO Ksyn Off: display bug fixed.
- Finish MIDI USB (MIDI out).
- Include other filter from preenfm2 forum.
- Add filter to perf page.
- Bug fix: ALGO 8, Mix/Pan5 in algo was modified by Mix/Pan2 on the preenfm2.
- Improve Rand menu.
- Volume (voice/timbre/out) is now OK. Oscillo OK per timbre.
- Put back play note from keyboard.
- Remove 2 of the 6 perfs.
- Go back to mixer (not edit) after loading mixer or preset.
- Bootloader.
- Save combo/mixer.
- Put back Scala per instrument.
- Init Scala scales after mixer load.
- Scala per instrument.
- Menu: move back from regular button to isolated one.
- Get rid of allChars. Use direct chars.
- Use mixer range info (MIDI first, last, shift).
- Use number of voices of mixer (prevent inconsistency with poly/mono).
- Instrument engine: replace number of voices per poly/mono.
- MIDI channel: add "All" option.
- Nicer (slightly bigger) buttons. More compact chars.
- MIDI general settings + global tuning.
- Display refactoring (mixer, editor, menu).
- Create smaller font to display MIDI activity and MIDI clock.
- Get rid of SysEx?
- UI: show currently edited OP in the algo.
- UI: show IM modified in the algo (1 sec).
- EPIC: menu button.
- Add all other algos.
- Remove SynthParamChecker.h.
- Colors: mixer in green, editor in blue, menu in red.
- Editor => OP (Osc): values are not correct.
- Make editor value yellow when changing (like in mixer).
- Encoder: add long press (0.5 sec).
- Editor: long press goes back to main editor page.
- Display note for MIDI clock.
- Bug: modified preset "*" wrong background.

