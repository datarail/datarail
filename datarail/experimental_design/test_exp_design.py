from process_assay import read_input
from designer import make_layout

csv_file = 'test_user_input2.csv'
plate_dims = [16, 24]
barcode_prefix = 'DUBI1_'
num_replicates = 3

treatments_df = read_input(csv_file, plate_dims,
                           barcode_prefix, encode_plate=True,
                           num_replicates=3)

drug_treatments = treatments_df[0]
nc_treatments = treatments_df[1]
pc_treatments = treatments_df[2]
bc_treatments = treatments_df[3]


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

Designs = make_layout(treatments_df, barcode_prefix,
                      encode_barcode=True,
                      plate_dims=[16, 24], nreps=num_replicates,
                      randomize=True, edge_bias=True)
