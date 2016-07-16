from datarail.utils.plate_fcts import *
import pandas as pd
import numpy as np

# define the source colors(drugs)
colors = ['Blue', 'Red', 'Green']

# define the destination plate
keys = ('well', 'drug', 'concentration')
dest_plt = pd.DataFrame([['A01', colors[0], 1]],
                                columns=keys)


# letter T
for i in range(0,7):
    dest_plt = dest_plt.append(pd.DataFrame([[chr(67+i)+'08', colors[i%3], 1]],
                                columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame([[chr(66)+'%02i'%(i+5), colors[i%3], 1]],
                                             columns=keys))

# letter L
for i in range(0,7):
    dest_plt = dest_plt.append(pd.DataFrame([[chr(66+i)+'13', colors[0], 1]],
                                columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame([[chr(73)+'%02i'%(i+13), colors[1], 1]],
                                             columns=keys))

# letter O
for i in range(0,5):
    dest_plt = dest_plt.append(pd.DataFrame([[chr(75+i)+'02', colors[1], 1]],
                                columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame([[chr(74)+'%02i'%(i+2), colors[2], 1]],
                                             columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame([[chr(75+i)+'06', colors[1], 1]],
                                columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame([[chr(80)+'%02i'%(i+2), colors[2], 1]],
                                             columns=keys))


# letter X
for i in range(0,7):
    dest_plt = dest_plt.append(pd.DataFrame([[chr(74+i)+'%02i'%(13+i), colors[i%3], 1]],
                                columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame([[chr(82-i)+'%02i'%(i+13), colors[i%3], 1]],
                                             columns=keys))



nmax = 4. # a new well in the source will be created every four dispense
col_cnt = [0]*len(colors)
# define the source and destination plates
src_plt = pd.DataFrame([ [chr(ord('A')+i)+'01', c, 1] for i,c in enumerate(colors)], columns=keys)
for i in range(0,len(dest_plt)):


    idx = [j for j, c in enumerate(colors) if dest_plt.drug.iloc[i] is c][0]


    if col_cnt[idx]%nmax==0 and col_cnt[idx]>0:
        # add a new well in the source if needed
        src_plt = src_plt.append(pd.DataFrame([[chr(ord('A')+idx)+'%02i'%(1+col_cnt[idx]/nmax),
                                                colors[idx], 1+col_cnt[idx]/nmax]],
                                              columns=keys))

    col_cnt[idx] += 1

    dest_plt.concentration.iloc[i] = np.ceil(col_cnt[idx]/nmax)
