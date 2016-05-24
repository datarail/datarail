# DATARAIL experimental design

## Class/structure for reagents, plates, ...

### Drug class
Contains information about a drug:
* Name
* Display name
* database ID (including batch)
* stock concentration

### Perturbation class
Generic class for information about experimental perturbations which are not a treatment
such as:
* cell line plated
* seeding density
* well volume
* 

### Plate class
Contains the experimental design for a plate (xarray):
* Design Name
* plate size
* treated wells
* matrix with vehicle
* treatments (multiple matrices with corresponding reagent)


## Function to design plate layout and store structures using python xarray

### Layout single agents
Function to layout single agent dose responses. Inputs are:
* list of drugs (structure format)
* list of concentrations

### Layout combination treatments
Function to layout cartesian combinations of drug treatmemts. Inputs are:
* list of drugs (structure format)
* list of concentrations
* list of drug pairs to test

  