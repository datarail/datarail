import functions as fct
import pandas as pd
import numpy as np
reload(fct)

# create the plate info file

df_plates = pd.DataFrame(dict({
    'barcode':['MH1_%02i'%i for i in range(1,21)],
    'cell_line':['CL_%1i'%(i%5+1) for i in range(0,20)],
    'treatment_file':['-']*5+['Example1_treatment_%1i.tsv'%i for i in [1]*5+[2]*5+[3]*5],
    'treatment_duration':[0]*5+[72]*15}))


df_plates.to_csv('./OUTPUT/Example1_plate_info.tsv', index=False, sep='\t')


# create the treatment files
np.random.seed(0)

df_trtfiles = [];

for i in range(3):
    well, agent, concentration, role = [[] for _ in range(4)]

    cnt = 0
    for ir in range(3,15):
        for ic in range(3,23):
            well += [chr(64+ir)+str(ic).zfill(2)]
            if cnt<20:
                agent += ['-']
                concentration += [0]
                role += ['negative_control']
            else:
                agent += ['D_' + str(np.mod(cnt/9,8)+1)]
                concentration += [10**(-3+.5*np.mod(cnt,9))]
                role += ['treatment']
            cnt += 1

    df_trtfiles += [pd.DataFrame(dict({
        'well':[well[j] for j in np.random.permutation(len(well))],
        'agent':agent,
        'concentration':concentration,
        'role':role})).sort_values('well')]

    df_trtfiles[i].to_csv('./OUTPUT/Example1_treatment_' + str(i+1) + '.tsv', index=False, sep='\t')


# create the output files with cell counts

df_treatments = pd.DataFrame([])
for i  in range(0,len(df_plates)):
    if df_plates.treatment_file[i] is '-':
        # no treatment
        df = df_trtfiles[0].copy()
        df['agent'] = '-'
        df['concentration'] = 0
        df['role'] = 'untreated'
    else:
        trtidx = int(df_plates.treatment_file[i][19])
        df = df_trtfiles[trtidx-1].copy()
    df['barcode'] = df_plates.barcode[i]
    df['cell_line'] = df_plates.cell_line[i]
    df['treatment_duration'] = df_plates.treatment_duration[i]
    df_treatments = df_treatments.append(df)

df_treatments.reset_index(inplace=True, drop=True)
df_sensitivity = pd.read_csv('./INPUT/agent_response.tsv', delimiter='\t')
df_growth = pd.read_csv('./INPUT/cell_line__growth.tsv', delimiter='\t')

df_res = fct.multiple_responses(df_treatments, df_growth, df_sensitivity)
df_bias = fct.low_cell_count_plate(fct.add_edge_effect(df_res, .07),
                                   ['MH1_10', 'MH1_18', 'MH1_19'], .9)

df_Col = fct.Columbus_fixed(df_treatments, df_bias)

df_res.to_csv('./OUTPUT/Example1_unbiased_results.tsv', index=False, sep='\t')
df_bias.to_csv('./OUTPUT/Example1_biased_results.tsv', index=False, sep='\t')
df_Col.to_csv('./OUTPUT/Example1_Columbus_output.tsv', index=False, sep='\t')
