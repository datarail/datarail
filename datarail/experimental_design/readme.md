# DATARAIL experimental design

Scripts for designing an experiment and recording metadata. This is not covering processing results or merging experimental design with data.

### Functions for plate layouts
High-level functions to design plate layouts. The treatments are constructed iteratively by adding conditions in a pandas table which is then laid on a plate (plate class).

#### Add conditions to treatment table with cartesian combinations of drugs and their concentrations
Inputs are:
* list of drugs (n classes)
* list of concentrations (list of 1 x n arrays)
* ?flag to define the order of combinations (or the highest order)?
This can be used for single agent dose-response, pairwise, and higher order combinations 

#### Add a user-defined condition to treatment table 
Input is:
* list of tuples of drug class and concentration

#### Lay out treatments on plate
Inputs are:
* flag on how to handle the controls (random, edge layout)
* seed for randomization
* number of randomized layouts
 
#### Add a reagent on specific wells
Inputs are:
* reagent (class)
* concentration layout (2D array)

#### Add a perturbation on specific wells
Inputs are:
* reagent (class)
* concentration layout (2D array)

### Functions for designing experiments and saving records
High-level functions to desing experiments such as define treated plates, combine treatment plates and treated plates.

### Additional functions

#### Create instances of `Reagent` class

#### Fetch reagent information from database
Inputs are:
* name or identifier
* url/link to api, or flag for internal database


### Classes
Classes used in the experimental design.

#### Reagent class
Contains general information about a reagent (see also subclasses):
* Display name (for user)
* Name (as in reference database)
* Database ID (including batch ID)
* Source database (with url)
* Comment (free text)

##### Drug subclass
Based on the class `Reagent`
* Stock concentration
* Vehicle (DMSO, aqueous, ...)

##### Cell line subclass
Based on the class `Reagent`
* Passage number

#### Perturbation class
Generic class for information about experimental perturbations which are not a treatment
such as cell line plated, seeding density, well volume.
* Name
* Layout (2D array with boolean, strings, numeric values, or reagent classes)
* Comment

#### Plate design class
Contains the experimental design for a plate (xarray):
* Design Name
* Plate size
* Treated wells
* Matrix with vehicle
* Treatments (numerical matrices with information from reagent class)
* Perturbations (2D array with boolean, strings, numeric values, or reagent classes. The layer has a name and can have a comment)

#### Treated plate class
Contains the metadata for a plate used in an experiment:
* Barcode
* Timestamps for plating, treatment, fixing, ...

#### Experiment class
Contains an instance of an experiment (metadata, design, timestamp):
* Mapping from treated plates to plate designs
* Comments



  
