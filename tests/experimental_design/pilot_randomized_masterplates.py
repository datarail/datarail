from datarail.utils.plate_fcts import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm


# define the source colors(drugs)
colors = ['Blue', 'Red', 'Green']

# define the destination plate
keys = ('well', 'drug', 'concentration')
dest_plt = pd.DataFrame([['A01', colors[0], 1]], columns=keys)


# letter T
for i in range(7):
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(67+i)+'08', colors[i % 3], 1]],
        columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(66)+'%02i' % (i+5), colors[i % 3], 1]],
        columns=keys))

# letter L
for i in range(7):
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(66+i)+'13', colors[0], 1]],
        columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(73)+'%02i' % (i+13), colors[1], 1]],
        columns=keys))

# letter O
for i in range(5):
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(75+i)+'02', colors[1], 1]],
        columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(74)+'%02i' % (i+2), colors[2], 1]],
        columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(75+i)+'06', colors[1], 1]],
        columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(80)+'%02i' % (i+2), colors[2], 1]],
        columns=keys))


# letter X
for i in range(7):
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(74+i)+'%02i' % (13+i), colors[i % 3], 1]],
        columns=keys))
    dest_plt = dest_plt.append(pd.DataFrame(
        [[chr(80-i)+'%02i' % (i+13), colors[i % 3], 1]],
        columns=keys))

nmax = 4.  # a new well in the source will be created every four dispense
col_cnt = [0]*len(colors)
# define the source and destination plates
src_plt = pd.DataFrame([[chr(ord('A')+i)+'01', c, 1]
                        for i, c in enumerate(colors)], columns=keys)
for i in range(len(dest_plt)):
    idx = [j for j, c in enumerate(colors)
           if dest_plt.drug.iloc[i] is c][0]

    if col_cnt[idx] % nmax == 0 and col_cnt[idx] > 0:
        # add a new well in the source if needed
        src_plt = src_plt.append(pd.DataFrame(
            [[chr(ord('A')+idx)+'%02i' % (1+col_cnt[idx]/nmax),
              colors[idx], 1+col_cnt[idx]/nmax]], columns=keys))

    col_cnt[idx] += 1

    dest_plt.concentration.iloc[i] = np.ceil(col_cnt[idx]/nmax)

src_well = []
dest_well = []

concentrations = range(1, 8)
for color in colors:
    for c in concentrations:
        s_well = src_plt.well[(src_plt.drug == color) & (
            src_plt.concentration == c)].values
        if len(s_well) > 0:
            d_well = dest_plt.well[(dest_plt.drug == color) & (
                dest_plt.concentration == c)].values.tolist()
            src_well += [s_well[0]]*len(d_well)
            dest_well += d_well
df = pd.DataFrame(zip(src_well, dest_well),
                  columns=['Source_well', 'Destination_well'])

df.to_csv('OUTPUT/pilot_randomized_mapping.csv', index=False)


def split_text(s):
    from itertools import groupby
    for k, g in groupby(s, str.isalpha):
        yield ''.join(list(g))

arr = np.zeros([16, 24])

dest_plt.index = range(63)
for i in dest_plt.index:
    r, c = split_text(dest_plt.well.ix[i])
    r = ord(r) - 65
    c = int(c) - 1
    if dest_plt.drug.ix[i] == 'Red':
        arr[r, c] = 1
    elif dest_plt.drug.ix[i] == 'Green':
        arr[r, c] = 2
    elif dest_plt.drug.ix[i] == 'Blue':
        arr[r, c] = 3

cmap = ListedColormap(['white', 'red', 'blue', 'green'])
bounds = [0, 1, 2, 3, 4]
norm = BoundaryNorm(bounds, cmap.N)
plt.pcolor(arr, cmap=cmap, norm=norm, edgecolor='k')
plt.gca().invert_yaxis()
plt.xticks(np.arange(0, 24)+0.5, [i for i in range(1, 25)])
plt.xlim([0, 24])
plt.yticks([i + 0.5 for i in range(16)],
           [chr(65+i) for i in range(16)])
plt.show()
