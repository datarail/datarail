import sys
sys.path.append('../..')
from datarail.databases import lincs_client

class Drug(object):
    def __init__(self, name, stock_conc=None):
        self.name = name
        # self.database_name = lincs_client.get_name(name)
        self.database_id = lincs_client.get_drug_id([name])
        self.url = 'http://lincs.hms.harvard.edu/db/sm/?search=&output_type=.csv'

        if stock_conc is None:
            stock_conc = 'User has not provided stock concentration'

        self.stock_concentration = stock_conc


    def __repr__(self):
        return self.name


class Cell_line(object):
    def __init__(self, name, vehicle=None, seeding_density=None):
        self.name = name
        self.database_id = lincs_client.get_cell_line_id([name])
        self.url = 'http://lincs.hms.harvard.edu/db/cells/?search=&output_type=.csv'

        if vehicle is None:
            vehicle = 'User has not provided vehicle information'
        if seeding_density is None:
            seeding_density = 'User has not provided seeding density'

        self.vehicle = vehicle
        self.seeding_density = seeding_density


    def __repr__(self):
        return self.name


    # def stock_concentration(self, stock_conc):
    #     self.stock_concentraion = stock_conc

        
    # def seeding_density(self, reagent_user_input):
    #     self.seeding_density = seeding_density


    # def vehicle(self, vehicle):
    #     self.vehicle = vehicle
        

# class Perturbation(object):
#     def __init__(self, perturbation_user_input):
#         self.name =
#         self.layout =  "2D layout "
#         self.comment = "TBD"


# class Plate_design(object):
#     '''will contain the Xaray structure ''' 
#     def __init__(self):
#         self.design_name
#         self.plate_size
#         self.treated_wells
#         self.matrix_vehicle
#         self.treatments =
#         self.pertubrations =
        

#  class treated_plate(object):
#      ''' contains barcode and timestamps for plating, treatment '''

#  class Experiment_class(object):
#      ''' will contan mapping between design and actual experiment ''''

 
