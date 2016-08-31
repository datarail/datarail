import pandas as pd
from collections import OrderedDict
from well_mapper import get_well_name


def to_tsv(Designs):

    plates = Designs.coords['plates'].values.tolist()

    for plate in plates:
        file_name = "%s_design.tsv" % plate
        design = Designs.sel(plates=plate)
        df = make_tsv(design, file_name)
    return df


def make_tsv(design, file_name):

    dimension_fields = design.dims.keys()

    coordinate_fields = [cf for cf in design.coords.keys()
                         if cf != 'plates']

    metadata_fields = [mf for mf in coordinate_fields
                       if mf not in dimension_fields]

    data_variables = design.data_vars.keys()

    well_info = OrderedDict()

    plate_dims = [16, 24]
    plate_size = 384
    well_info['well'] = [get_well_name(w, plate_dims)
                         for w in range(plate_size)]

    for field in metadata_fields:
        field_array = design[field].values.reshape(1, 384)[0]
        well_info[field] = field_array.tolist()

    for field in data_variables:
        field_array = design[field].values.reshape(1, 384)[0]
        well_info[field] = field_array.tolist()

    df = pd.DataFrame(well_info)
    df.to_csv(file_name, sep='\t')
    return df
