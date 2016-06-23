from functions import *


# create the plate info file

df_plates = pd.DataFrame(dict({
    'barcode':['20160623_MH_%02i'%i for i in range(1,13)],
    'cell_line':['CL_%1i'%(i%6+1) for i in range(0,12)],
    'treatment_file':['ExampleTrt_%1i.tsv'%i for i in [1]*6+[2]*6],
    'treatment_duration':[72]*12}))


df_plates.to_csv('../OUTPUT/Example1_plate_info.tsv', index=False, sep='\t')


# create the treatment files
np.random.seed(0)

df_trtfiles = [];

for i in (1,2):
    well = []
    agent = []
    concentration = []

    for ir in range(2,16):
        for ic in range(3,23):
            well += [chr(63+ir)+str(ic).zfill(2)]
            if np.random.rand()<.1:
                agent += ['-']
                concentration += [0]
            else:
                agent += ['D_' + str(np.random.randint(4)+1)]
                concentration += [10**(-3+.5*np.random.randint(9))]

    df_trtfiles += [pd.DataFrame(dict({
        'well':well,
        'agent':agent,
        'concentration':concentration}))]

    df_trtfiles[i-1].to_csv('../OUTPUT/Example1_treatment' + str(i) + '.tsv', index=False, sep='\t')


# creat the output files with cell counts

df_treatments = pd.DataFrame([])
for i  in range(0,len(df_plates)):
    trtidx = int(df_plates.treatment_file[2][11])
    df = df_trtfiles[trtidx-1].copy()
    df['barcode'] = df_plates.barcode[i]
    df['cell_line'] = df_plates.cell_line[i]
    df['treatment_duration'] = df_plates.treatment_duration[i]
    df_treatments = df_treatments.append(df)

df_treatments.reset_index(inplace=True)
df_sensitivity = pd.read_csv('../INPUT/agent_response.tsv', delimiter='\t')
df_growth = pd.read_csv('../INPUT/cell_line__growth.tsv', delimiter='\t')

df_res = multiple_responses(df_treatments, df_growth, df_sensitivity)
df_bias = add_edge_effect(df_res, .1)

df_Col = Columbus_fixed(df_treatments, df_bias)

df_res.to_csv('../OUTPUT/Example1_unbiased_results.tsv', index=False, sep='\t')
df_bias.to_csv('../OUTPUT/Example1_biased_results.tsv', index=False, sep='\t')
df_Col.to_csv('../OUTPUT/Example1_Columbus_output.tsv', index=False, sep='\t')
