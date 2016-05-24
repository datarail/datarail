# DATARAIL experimental design

Scripts for designing an experiment and recording metadata. This is not covering processing results or merging experimental design with data.

### Functions for plate layouts
High-level functions to design plate layouts.

#### Layout single agents dose responses
Inputs are:
* list of drugs (n classes)
* list of concentrations (list of 1 x n arrays)

#### Layout cartesian combinations of drug treatmemts
Inputs are:
* list of drugs (n classes)
* list of concentrations (list of 1 x n arrays)
* list of drug pairs to test (2 x m array with values in 1..n)
 
#### Add a reagent on specific wells
Inputs are:
* reagent (class)
* concentration layout (2D array)

#### Randomize positions with defined control positions
Inputs are:
* flag on how to handle the controls (random, edge layout)

### Functions for desinging experiments and saving records
High-level functions to desing experiments.


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
* Seeding density

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
* Treatments (multiple matrices with corresponding reagent)
* Perturbations

#### Treated plate class
Contains the metadata for a plate used in an experiment:
* Barcode
* Timestamps for plating, treatment, fixing, ...

#### Experiment class
Contains an instance of an experiment (metadata, design, timestamp):
* Mapping from treated plates to plate designs
* Comments



  
