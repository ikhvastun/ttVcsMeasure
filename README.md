# ttVcsMeasure
The code is used to run the ttV cross section measurement with 2016 and 2017 datasets

## Compilation
To compile the code use [MAKE](http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/)

```make
make -f makefiles/makeTTZ clean; make -f makefiles/makeTTZ
```
## Running 2nd skimming
Before producing nice plots, the one should run 2nd skimming on the tuples produced by the [ntuplizer](https://github.com/GhentAnalysis/heavyNeutrino). In this 2nd skimming, the one can put tighter selection on the leptons, jets, invariant mass of 2 leptons, etc. To run 2nd skimming:

```bash
./skimTree [one file to run on] 
```

Use `skimTuples.sh` from `scripts/` folder to submit on the whole sample

```bash
bash skimTuples.sh 
```

The output of the latter command is a bunch of small root files, which the one should combine in an aggregate root file for a particular sample

## Plots/Tables/Datacards production for particular analysis 
To run the code to produce plots/tables/datacards for ttZ cross section measurement, use the following semantics: 

```bash
./main [files to run on] [options to run] [selections you wanna use] [process to run on] [event number]
```

The following possible options are:

- [files to run on]: use 1 or 2 files from `data/samples/` directory, if you wanna use 2 files use `,` as a separator
```bash
./main data/samples/samples_ttZ_npDD.txt,data/samples/samples_ttZ_2017_npDD.txt runFullSelection selection:ttZclean
```

- [options to run]: (runFullSelection, runOnOneProcess, debug)
  * runFullSelection -> runs the whole code over all files specified in previous option
  * runOnOneProcess -> run on one specific process, for example ttZ, ttX, etc. If you wanna use it, put the process to run on in the end of the command

  ```bash
  ./main data/samples/samples_ttZ_npDD.txt runOnOneProcess selection:ttZclean ttZ
  ```

  * debug -> shows all possible information (lepton pt, jet pt, etc), in this case specify in the end process to run on and event number

  ```bash
  ./main data/samples/samples_ttZ_npDD.txt debug selection:ttZ data 1594370560
  ```

- [selections you wanna use] : specify selection (you wanna look into control region WZ or ttZ enriched region). Put selection: before the name of the selection. Possible options are: (ttZ3L, ttZ4L, ttZ, ttZclean, WZ, ZZ, ttbar, DY, Xgamma)

To submit code on T2 use `submitCR.sh` code in `scripts/` directory.

```bash
sh submitCR.sh
```


