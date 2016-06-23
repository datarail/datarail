import numpy as np

def trimmean(vals):
    return np.mean(sorted(vals)[int(len(vals)*.25):int(len(vals)*.75)])
