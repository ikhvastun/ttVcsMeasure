import numpy as np
fnameM = '../data/dumps/events_data_3mu_2017_test.txt'
#fnameW = '/user/wverbeke/dump.txt'
fnameW = 'dumps/daniele_3mu_data_2017.txt'
eventM = np.loadtxt(fnameM, dtype = int)
#eventW = np.loadtxt(fnameW, dtype = int, delimiter = " ", usecols = (0,1,2))
eventW = np.loadtxt(fnameW, dtype = int)

eventM_hashable = map(tuple, eventM)
eventM = set(eventM_hashable)

eventW_hashable = map(tuple, eventW)
eventW = set(eventW_hashable)

print("events Daniel selects and I don't: " + str(len(eventW - eventM)))
print(eventW - eventM)
print("events I select and Daniel doesn't: " + str(len(eventM - eventW)))
print(eventM - eventW)
print("events in common: " + str(len(eventM & eventW)))