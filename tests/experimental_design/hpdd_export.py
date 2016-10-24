import sys
import os.path as op
import numpy as np
import pandas as pd
from datarail.experimental_design.hpdd_utils import export_hpdd
from datarail.experimental_design.well_mapper import get_well_name

design = pd.DataFrame()
design['plate'] = ['plate_1'] * 96 + ['plate_2'] * 96
design['well'] = [get_well_name(i, (8, 12)) for i in range(96)] * 2

agents = np.array([['Fluid %d' % i] * 12 for i in range(1, 8+1)]).flatten()
design['agent'] = np.tile(agents, 2)
concs = np.logspace(1, -2, 12)
design['concentration'] = np.tile(concs, 8 * 2)
dmso_control_loc = (design.plate == 'plate_1') & (design.well == 'H12')
design.loc[dmso_control_loc, 'agent'] = np.nan
design.loc[dmso_control_loc, 'concentration'] = np.nan

design['agent_2'] = np.nan
design['concentration_2'] = np.nan
f9_loc = (design.plate == 'plate_2') & (design.well.str.contains(r'01$'))
design.loc[f9_loc, 'agent_2'] = 'Fluid 9'
f9_concs = np.logspace(np.log10(7), np.log10(.007), 8)
design.loc[f9_loc, 'concentration_2'] = f9_concs

agents = pd.DataFrame()
agents['name'] = ['Fluid %d' % i for i in range(1, 9 +1)]
agents['stock_concentration'] = 10000

output_location = op.join(op.abspath(op.dirname(__file__)), 'OUTPUT')
output_path = op.join(output_location, 'test_example.hpdd')

export_hpdd(design, agents, output_path, assay_volume=100)

print "Output saved to:", output_path
