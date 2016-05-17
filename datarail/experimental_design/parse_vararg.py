def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='parse input feild for Deisgn plate')

    parser.add_argument('-ComboDoses',
                        default=[],
                        type=lambda x:
                        all([isinstance(y, (int, float)) for y in x]))
    parser.add_argument('-ComboLists',
                        default=[],
                        type=lambda x:
                        all([isinstance(y, (int, float)) for y in x]))
    parser.add_argument('-DrugPairs',
                        default=0,
                        type=lambda x: isinstance(x, (int, float)))
    parser.add_argument('-Seed', default=1,
                        type=lambda x: isinstance(x, (int, float)))
    parser.add_argument('-edge_ctrl',
                        default=True,
                        type=lambda x: isinstance(x, bool))
    parser.add_argument('-stock_conc',
                        default=1e4,
                        type=lambda x: isinstance(x, (int, float)))
    parser.add_argument('-Vehicle',
                        default='DMSO')
    parser.add_argument('-well_volume',
                        default=60,
                        type=lambda x: isinstance(x, (int, float)) &
                        (not isinstance(x, list)))
    parser.add_argument('-plate_dims',
                        default=[16, 24],
                        type=lambda x:
                        all([isinstance(y, (int, float)) for y in x]) &
                        (isinstance(x, list)))
    parser.add_argument('-Perturbations')

    # Based on the specifications of the D300
    parser.add_argument('-min_volume',
                        default=1.3e-5,
                        type=lambda x: isinstance(x, (int, float)))
    parser.add_argument('-step_volume',
                        default=2e-5,
                        type=lambda x: isinstance(x, (int, float)))
    parser.add_argument('-max_DMSOpct', default=0.2,
                        type=lambda x: isinstance(x, (int, float)))

    args = parser.parse_args()
    return args

   
