import pandas as pd
import numpy as np
import datetime


######  functions to generate synthetic drug response data  #######

# sigmoidal function for generating cell count (concentrations in log10 domain)
def cellcount(conc, max_count, EC50, Hill, min_count):
        return min_count + (max_count - min_count)/(1. + ((10.**conc)/(10.**EC50))**Hill)


# sigmoidal function for generating dead cell count (concentrations in log10 domain)
def deadcount(conc, max_count, EC50, Hill, min_count):
        return max_count - (max_count - min_count)/(1. + ((10.**conc)/(10.**EC50))**Hill)

# add binomial noise 
def add_count_noise(cc, coeff=0.5):
    if type(cc) is list or float:
        cc = np.asarray(cc)
        
    return np.ceil(np.random.binomial((cc/coeff).astype(int), coeff))

# generate drug response data value
def allcounts(conc, x0, xctrl, GRmax, x50, Hill):
    xend = x0 * 2**(np.log2(GRmax+1.)*np.log2(1.*xctrl/x0))
    dead_max = (xend+x0/10)*max(0., min(0.2-GRmax,.9))/(1-max(0., min(0.2-GRmax,.6)))
    dead_min = np.ceil(x0/20.)

    totalcells = add_count_noise(cellcount(conc, xctrl+dead_min/2, x50-.1, Hill, xend+dead_max/2))
    deadcells = add_count_noise(deadcount(conc, dead_max, x50+.1, Hill, dead_min))/2 # dead cells are split with corpses
    viablecells = totalcells - deadcells
    corpses = add_count_noise(np.maximum(10, x0-totalcells)+10)+deadcells

    
    return pd.DataFrame([[conc, totalcells, deadcells, viablecells, corpses, (corpses+deadcells)/(corpses+totalcells)]],
                        columns=['concentration', 'cell_count__total','cell_count__dead','cell_count','corpse_count','frac_dead']) 
       # return pd.DataFrame(dict=({conc, totalcells, deadcells, viablecells, corpses],
        #                columns=['concentration', 'cell_count__total','cell_count__dead','cell_count','corpse_count']) 


# full dose response curve
def dose_response(x0, xctrl, GRmax, x50, Hill):
    conc = np.array(range(-7,4))/2.
    return allcounts(conc, x0, xctrl, GRmax, x50, Hill)


# values based on predefined table and selected conditions
def single_condition_response(df_sensitivity, df_growth, cell_line,  agent, concentration, time=72, x0=500):

    para = df_sensitivity[(df_sensitivity.agent==agent) & (df_sensitivity.cell_line==cell_line)]
    xctrl = np.ceil(x0 * 2**(time/df_growth[df_growth.cell_line==cell_line]['doubling_time'].values[0]))

    return allcounts(concentration, x0, xctrl, para['MaxEffect'].values[0],
                     para['log10(EC50)'].values[0], para['Hill'].values[0])
    


######### functions to generate artefacts #############

# edge effects
#def add_edge_effect(df_):
        


#########  reproduce outputs from microscopes  #########

# Columbus data
def Columbus_row_data(barcode, timestamp, well, df_values, Nfields=6):
    df_pos = pd.DataFrame([[barcode + ' > ' + timestamp,
                            well_as_1digit(well),
                            well_as_row_col(well_as_1digit(well))['row'],
                            well_as_row_col(well_as_1digit(well))['column']]],
                            columns=['Result', 'Well', 'Row', 'Column'])
    df_url = pd.DataFrame([[Nfields,
                            'http://columbus.hms.harvard.edu/browse/measurement/16514/well='+\
                            df_pos['Row'][0]+'.'+df_pos['Column'][0]]],
                          columns=['Number of Analyzed Fields', 'URL'])
    return pd.concat([df_pos, df_values, df_url], axis=1) 
    

def Columbus_fixed(df_trt, df_gr, df_sens):
    # usage: Columbus_fixed(ref['df_trt'], ref['df_gr'], ref['df_sens'])
    df_res = pd.DataFrame([])
    timestamp = '2016-06-06 12:34:56'

    for i in range(0,len(df_trt)):
        res = single_condition_response(df_sens, df_gr, df_trt['cell_line'][i],
                    df_trt['agent'][i], df_trt['concentration'][i])
        df_ = Columbus_row_data(df_trt['barcode'][i], timestamp,
                                        df_trt['well'][i], res, 5+1.*(np.random.uniform()>.1))
        df_res = df_res.append(df_)

    df_res.drop(['concentration', 'cell_count', 'frac_dead'], 1, inplace=True)
    df_res.rename(index=str, columns={'cell_count__total': 'Hoechst_pos',
                                      'cell_count__dead': 'Hoechst_LDR_pos',
                                      'corpse_count': 'LDR_pos_Hoechst_neg'},
                  inplace=True)
    return df_res

def Columbus_timecourse(df_trt, df_gr, df_sens):
    # usage: Columbus_fixed(ref['df_trt'], ref['df_gr'], ref['df_sens'])
    df_res = pd.DataFrame([])

    barcodes = df_trt.barcode.drop_duplicates()
    for it in range(0, 84, 4):        
        for i in range(0,len(df_trt)):
            barcode_idx = [j for j in range(0,len(barcodes)) if df_trt.barcode[i] is barcodes.iloc[j]][0]
            timestamp = (datetime.datetime(2016, 6, 6, 9, 10, 11) +\
                         datetime.timedelta(hours=it, minutes=20*barcode_idx)).strftime('%Y-%m-%d %H:%M:%S')

            res = single_condition_response(df_sens, df_gr, df_trt['cell_line'][i],
                                            df_trt['agent'][i], df_trt['concentration'][i])
            df_ = Columbus_row_data(df_trt['barcode'][i], timestamp,
                                    df_trt['well'][i], res, 5+1.*(np.random.uniform()>.1))
            df_res = df_res.append(df_)

    df_res.drop(['concentration', 'cell_count'], 1, inplace=True)
    df_res.rename(index=str, columns={'cell_count__total': 'Hoechst pos - LDR neg',
                                      'cell_count__dead': 'Hoechst pos - LDR pos',
                                      'cell_count__corpse': 'Hoechst neg - LDR pos'},
                  inplace=True)
    return df_res

######## initialization ########
def load_references():
    df_growth = pd.read_csv('../INPUT/cell_line__growth.tsv', delimiter='\t')
    df_sensitivity = pd.read_csv('../INPUT/agent_response.tsv', delimiter='\t')
    df_treatments = pd.read_csv('../INPUT/treatment_list.tsv', delimiter='\t')
    return dict({'df_gr':df_growth, 'df_sens':df_sensitivity, 'df_trt':df_treatments})
           


#########  misc functions ############

def trimmean(vals):
    return np.mean(sorted(vals)[int(len(vals)*.25):int(len(vals)*.75)])


def well_as_1digit(well):
    if well[1]=='0':
        well = well[0]+well[2]
    return well

def well_as_row_col(well):
    return dict({'row':chr(ord(well[0])-16),'column':well[1:]})

                    
