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

    
    return pd.DataFrame([[totalcells, deadcells, viablecells, corpses, (corpses+deadcells)/(corpses+totalcells)]],
                        columns=['cell_count__total','cell_count__dead','cell_count','corpse_count','frac_dead']) 
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

    if len(para)==0:
        return allcounts(concentration, x0, xctrl, 1, np.inf, .01)
    else:
        return allcounts(concentration, x0, xctrl, para['MaxEffect'].values[0],
                         para['log10(EC50)'].values[0], para['Hill'].values[0])
    

# generate all results based on the treatments, cell line division rates and drug sensitivity
# multiple_responses(ref['df_rtr'], ref['df_gr'], ref['df_sens'])
def multiple_responses(df_trt, df_gr, df_sens):
        
    df_res = pd.DataFrame([])

    for i in range(0,len(df_trt)):
        res = single_condition_response(df_sens, df_gr, df_trt['cell_line'][i],
                                        df_trt['agent'][i], df_trt['concentration'][i])
        res = pd.concat([df_trt.iloc[i:(i+1),:].reset_index(drop=True), res], axis=1)
        
        df_res = df_res.append(res)
        
    return df_res.reset_index(drop=True)


######### functions to generate artefacts #############

# edge effects :
#   - assume that plates have the format (2^n, 1.5*2^n)
#   - need variables 'Row' and 'Column'
#   - change some cells from 'viable' to 'dead'
def add_edge_effect(df_, strength=.3):

    added_Column = 'Column' not in df_.keys()
    if added_Column:
        df_ = pd.concat([df_, pd.DataFrame([well_as_row_col(df_.well.iloc[i])['column'] \
                                      for i in range(0,len(df_))], columns=['Column'])], axis=1)
        
    added_Row = 'Row' not in df_.keys()
    if added_Row:
        df_ = pd.concat([df_, pd.DataFrame([well_as_row_col(df_.well.iloc[i])['row'] \
                                      for i in range(0,len(df_))], columns=['Row'])], axis=1)
        
    Nrows = 2**(np.ceil(np.log2( max(df_.Row) )))
    Ncols = 1.5*2**(np.ceil(np.log2( max(df_.Column)/1.5 )))
    EdgeDist = np.minimum(
        np.minimum(Ncols-df_.Column,df_.Column),
        np.minimum(Nrows-df_.Row,df_.Row))

    # exponential decay, half-distance is 1.5
    Effect = strength*np.exp(-(EdgeDist-1)/1.5)
    NewDeadCells = np.ceil(df_.cell_count*Effect);
    df_.cell_count__dead += NewDeadCells
    df_.cell_count -= NewDeadCells
    df_.frac_dead = (df_.corpse_count+df_.cell_count__dead)/(df_.corpse_count+df_.cell_count__total)

    if added_Column:
        df_.drop('Column', 1, inplace=True)

    if added_Row:
        df_.drop('Row', 1, inplace=True)
        
    return df_


#########  reproduce outputs from microscopes  #########

defaultNfield = 6 # general scope variable

# generate a row data as exported from Columbus
def Columbus_row_data(barcode, timestamp, well, df_values, Nfields=6):
    df_pos = pd.DataFrame([[barcode + ' > ' + timestamp,
                            well_as_1digit(well),
                            well_as_row_col(well_as_1digit(well))['row'],
                            well_as_row_col(well_as_1digit(well))['column']]],
                            columns=['Result', 'Well', 'Row', 'Column'])
    df_url = pd.DataFrame([[Nfields,
                            'http://columbus.hms.harvard.edu/browse/measurement/16514/well='+\
                            str(df_pos['Row'][0])+'.'+str(df_pos['Column'][0])]],
                          columns=['Number of Analyzed Fields', 'URL'])
    
    # account for missing fields
    for f in ['cell_count__total', 'cell_count__dead', 'cell_count', 'corpse_count']:
        df_values[f] = (df_values[f]*Nfields)/defaultNfield

    if 'barcode' in df_values.columns:
        df_values.drop('barcode', 1, inplace=True)
    
    return pd.concat([df_pos, df_values, df_url], axis=1) 
    

# generate data in the export format of Columbus
# needs a list of treatments and the corresponding drug responses
# usage: Columbus_fixed(ref['df_trt'], df_res)
def Columbus_fixed(df_trt, df_res):
        
    timestamp = '2016-06-06 12:34:56'
    df_Columb = pd.DataFrame([])
    for i in range(0,len(df_trt)):
        df_ = Columbus_row_data(df_trt['barcode'][i], timestamp,
                                df_trt['well'][i], df_res.iloc[i:(i+1),:].reset_index(drop=True),
                                5+1.*(np.random.uniform()>.1))
        
        df_Columb = df_Columb.append(df_)
        
    df_Columb.drop(['concentration', 'cell_count', 'frac_dead', 'agent', 'cell_line', 'well'], 1, inplace=True)
    df_Columb.rename(index=str, columns={'cell_count__total': 'Hoechst_pos',
                                      'cell_count__dead': 'Hoechst_LDR_pos',
                                      'corpse_count': 'LDR_pos_Hoechst_neg'},
                  inplace=True)
    return df_Columb


# generate timecourse data in the export format of Columbus (t=0, 4, 8, ... 80h)
# needs a list of treatments, cell line division rates, and drug responses
# usage: Columbus_timecourse(ref['df_trt'], ref['df_gr'], ref['df_sens'])
def Columbus_timecourse(df_trt, df_gr, df_sens, strength=0):
    print 'check implementation --- not fully done --- '
    
    df_res = pd.DataFrame([])
    df_Columb = pd.DataFrame([])
    
    barcodes = df_trt.barcode.drop_duplicates()
    for it in range(0, 84, 4):        
        for i in range(0,len(df_trt)):
            barcode_idx = [j for j in range(0,len(barcodes)) if df_trt.barcode[i] is barcodes.iloc[j]][0]
            timestamp = (datetime.datetime(2016, 6, 6, 9, 10, 11) +\
                         datetime.timedelta(hours=it, minutes=20*barcode_idx)).strftime('%Y-%m-%d %H:%M:%S')
                        
            res = single_condition_response(df_sens, df_gr, df_trt['cell_line'][i],
                                            df_trt['agent'][i], df_trt['concentration'][i])
            df_res = df_res.append(res)
            
            df_ = Columbus_row_data(df_trt['barcode'][i], timestamp, df_trt['well'][i],
                                    add_edge_effect(res, strength),
                                    5+1.*(np.random.uniform()>.1))
            df_Columb = df_Columb.append(df_)

    df_res.drop(['concentration', 'cell_count'], 1, inplace=True)
    df_res.rename(index=str, columns={'cell_count__total': 'Hoechst pos - LDR neg',
                                      'cell_count__dead': 'Hoechst pos - LDR pos',
                                      'cell_count__corpse': 'Hoechst neg - LDR pos'},
                  inplace=True)
    return df_res



######## working example ########
# default example files
def load_references():
    df_growth = pd.read_csv('../INPUT/cell_line__growth.tsv', delimiter='\t')
    df_sensitivity = pd.read_csv('../INPUT/agent_response.tsv', delimiter='\t')
    df_treatments = pd.read_csv('../INPUT/MultiLines_perPlate_example.tsv', delimiter='\t')
    return dict({'df_gr':df_growth, 'df_sens':df_sensitivity, 'df_trt':df_treatments})
           

def Columbus_example(filename='example', strength=.3):
    ref = load_references()

    df_res = multiple_responses(ref['df_trt'], ref['df_gr'], ref['df_sens'])
    df_bias = add_edge_effect(df_res, strength)
    
    df_Col = Columbus_fixed(ref['df_trt'], df_bias)

    df_res.to_csv('../OUTPUT/unbiased_results_' + filename + '.tsv', index=False, sep='\t')
    if strength>0:
        df_bias.to_csv('../OUTPUT/biased_results_' + filename + '.tsv', index=False, sep='\t')
    df_Col.to_csv('../OUTPUT/Columbus_output_' + filename + '.tsv', index=False, sep='\t')
    
    return(df_res, df_bias, df_Col)

#########  misc functions ############

def trimmean(vals):
        
    return np.mean(sorted(vals)[int(len(vals)*.25):int(len(vals)*.75)])


def well_as_1digit(well):
    if well[1]=='0':
        well = well[0]+well[2]
    return well


def well_as_row_col(well):
    return dict({'row':ord(well[0])-64,'column':int(well[1:])})

                    
