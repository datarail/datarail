from process_assay import read_input
from designer import make_layout
import numpy as np

csv_file = 'test_user_input2.csv'
plate_dims = [16, 24]
barcode_prefix = 'DUBI1_'
num_replicates = 3

treatment_dicts = read_input(csv_file, plate_dims,
                             barcode_prefix, encode_plate=True,
                             num_replicates=3)


drug_treatments = treatment_dicts[0]
nc_treatments = treatment_dicts[1]
pc_treatments = treatment_dicts[2]
bc_treatments = treatment_dicts[3]

concentrations = np.hstack((1, 10 ** np.linspace(1, 4, 9)))
drug_class1 = ['IU1']
drug_class2 = ['MI-2']
combo_pairs = [(d1, d2) for d1 in drug_class1 for d2 in drug_class2]
combo_k = max([len(c) for c in combo_pairs])
combo_doses = {}
combo_drugs = drug_class1 + drug_class2
for drug in combo_drugs:
    combo_doses[drug] = concentrations[4:6]

print "\n Report of the experimental assay"
print "-----------------------------------\n"
print "This experiment assays %d drugs making up a total"\
    " of %d treatment conditions\n" % (len(drug_treatments.keys()),
                                       sum(len(v) for v
                                           in drug_treatments.itervalues()))


print "This experiment has %d negative control%s\n" % (
    len(nc_treatments.keys()), [p[0] if len(nc_treatments.keys()) ==
                                1 else p[1] for p in [('', 's')]][0])

print "This experiment has %d positive control%s\n" % (
    len(pc_treatments.keys()), [p[0] if len(pc_treatments.keys()) ==
                                1 else p[1] for p in [('', 's')]][0])

print "This experiment has %d compound%s used for encoding the barcode\n" % (
    len(bc_treatments.keys()), [p[0] if len(bc_treatments.keys()) ==
                                1 else p[1] for p in [('', 's')]][0])

Designs = make_layout(treatment_dicts, barcode_prefix,
                      encode_barcode=True,
                      plate_dims=[16, 24], nreps=num_replicates,
                      randomize=True, biased_randomization=True,
                      combo_pairs=combo_pairs, combo_doses=combo_doses,
                      combo_k=combo_k)
