# ttVcsMeasure
The code is used to run the ttV cross section measurement with 2016 and 2017 datasets

## Compilation
To compile the code use [MAKE](http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/)

```make
make -f makeTTZ clean; make -f makeTTZ
```
## Usage
To run the code use the following semantics: 

```bash
usage: git [--version] [--help] [-C <path>] [-c name=value]
./main [files to run on] [options to run] [selections you wanna use] [process to run on] [event number]
```

The following possible options are:

- [files to run on]: use 1 or 2 files from data/samples/ directory, if you wanna use 2 files use ',' as a separator
```bash
./main data/samples/samples_ttZ_npDD.txt,data/samples/samples_ttZ_2017_npDD.txt runFullSelection selection:ttZclean
```

- [options to run]: (runFullSelection, runOnOneProcess, debug)
   * runFullSelection -> runs the whole code over all files specified in previous option
   * runOnOneProcess -> run on one specific process, for example ttZ, ttX, etc. If you wanna use it, put the process to run on in the end of the command
```bash
./main data/samples/samples_ttZ_npDD.txt,data/samples/samples_ttZ_2017_npDD.txt runOnOneProcess selection:ttZclean ttZ
```
   * debug -> shows all possible information (lepton pt, jet pt, etc), in this case specify in the end process to run on and event number

```bash
./main data/samples/samples_ttZ_npDD.txt,data/samples/samples_ttZ_2017_npDD.txt debug selection:ttZ data 1594370560
```

- [selections you wanna use] : specify selection (you wanna look into control region WZ or ttZ enriched region). Put selection: before the name of the selection. Possible options are: (ttZ3L, ttZ4L, ttZ, ttZclean, WZ, ZZ, ttbar, DY, Xgamma)


