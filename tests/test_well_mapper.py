from datarail.experimental_design import well_mapper


def test_get_well_index():
    index = well_mapper.get_well_index('B10', [16, 24])
    assert index == 33


def test_get_well_name():
    name = well_mapper.get_well_name(33, [16, 24])
    assert name == 'B10'


    
