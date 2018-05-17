## Computational workflow for design of dose-response experiments

### Installation
* The repository can be installed from command line as shown below
   ```
   $ git clone https://github.com/datarail/datarail.git
   ```
* To install dependencies and enable importing modules from any location on your local machine, `cd` into the `datarail` folder, followed by the command below.
   ```
   $ pip install -e .
   ```

### Getting started
* Set up the well and plate level metadata files as shown in `datarail/examples`
* Start a Jupyter notebook or IPython session.
* The layout of drugs on doses across  96/384 well plates can be constructed using the code below. The pandas dataframe `dfm` contains the desingned layout. Refer to `datarail/examples` for a detailed explanation with examples. 
  ```python
  import pandas as pd
  from datarail.experimental_design import process_assay as pa
  dfp = pd.read_csv('plate_level_metadata.csv')
  dfm = pa.randomize_wells(dfp)
  dfm.to_csv('dose_response_layout_metadata.csv', index=False)
  ```
* The metadata file can be exported to a .hpdd file that can be used by the D300 printer. The stock concentraion for the drug also needs to be provided for each drug in the assay.
  ```python
  from datarail.experimental_design import hpdd_utils as hu
  hu.export_hpdd(dfm, dfs, 'dose_response_layout_metadata.hpdd')
  ```
* The layout can be visualized using the code below
  ```python
  from datarail.experimental_design import plot_plate_layout as ppl
  ppl.plot_summary(dfr, 'dose_response_layout_metadata.pdf')  
  ```
 * If the D300 software was used for designing the experiment, follow the steps below inorder to save the metadata file based on       `datarail` convention for subsequent downstream analysis.
   - Open the `.xml` from D300 in Excel and save as a `.xlsx` file.
   - Use the code below to save the metadata in a dataframe `dfm`
   ```python
   from datarail.experimental_design import export_D300_xml as edx
   dfm = edx.export2pd('D300_filename.xlsx')
   ```
  
