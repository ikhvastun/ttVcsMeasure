#!/bin/bash

# 2016
# 2L
#sh runCode.sh data/samples/samples_ttW_npDD_CMIDDD.txt ttW
#sh runCode.sh data/samples/samples_ttW_npDD_CMIDDD.txt ttWclean

# 3L
#sh runCode.sh data/samples/samples_ttZ_npDD.txt ttZ
#sh runCode.sh data/samples/samples_ttZ_npDD.txt ttZ3L
#sh runCode.sh data/samples/samples_ttZ_npDD.txt ttZ3Lclean

#sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_WZ.txt WZ
#sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_nonprompt.txt ttbar
#sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_nonprompt.txt DY
#sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_Zgamma.txt Zgamma

# 3L tZq
#sh runCode.sh data/samples/samples_ttZ_npDD.txt tZq

# 4L
#sh runCode.sh data/samples/samples_ttZ.txt ttZ4L
#sh runCode.sh data/samples/2016CR/samples_ttZ_ZZ.txt ZZ

###########################################################################
# 2017
# 2L
# in principle should be unblinded, but we look only in kinematic variables atm
#sh runCode.sh data/samples/samples_ttW_2017_npDD_CMIDDD.txt ttW

# 2017
# 3L
# not unblinded in 2017 yet
#sh runCode.sh data/samples/samples_ttZ_2017_npDD.txt ttZ3L

#sh runCode.sh data/samples/2017CR/samples_ttZ_2017_npDD_WZ.txt WZ
#sh runCode.sh data/samples/2017CR/samples_ttZ_2017_npDD_nonprompt.txt ttbar
#sh runCode.sh data/samples/2017CR/samples_ttZ_2017_npDD_nonprompt.txt DY
#sh runCode.sh data/samples/2017CR/samples_ttZ_2017_npDD_Zgamma.txt Zgamma

# 4L
# not unblinded in 2017 yet
#sh runCode.sh data/samples/samples_ttZ_2017.txt ttZ4L
#sh runCode.sh data/samples/2017CR/samples_ttZ_2017_ZZ.txt ZZ

###########################################################################
# combination of 2016 and 2017
# 3L
#sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_WZ.txt,data/samples/2017CR/samples_ttZ_2017_npDD_WZ.txt WZ
sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_nonprompt.txt,data/samples/2017CR/samples_ttZ_2017_npDD_nonprompt.txt ttbar
sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_nonprompt.txt,data/samples/2017CR/samples_ttZ_2017_npDD_nonprompt.txt DY
sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_Xgamma.txt,data/samples/2017CR/samples_ttZ_2017_npDD_Xgamma.txt Xgamma
sh runCode.sh data/samples/2016CR/samples_ttZ_npDD_ZZ.txt,data/samples/2017CR/samples_ttZ_2017_npDD_ZZ.txt ZZ
