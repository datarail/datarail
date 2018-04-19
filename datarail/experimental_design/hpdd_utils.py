"""Utilities for reading and writing Hewlett-Packard D300 .hpdd files."""

from __future__ import division
from lxml import objectify, etree
import numpy as np
from .. import __version__ as package_version
from .well_mapper import get_well_index


MU = u'\u00b5'


def _register_objectify_types():
    """Adjust XML serialization to match HPDD software output."""
    global _objectify_bool, _objectify_float
    # Serialize using upper case first letter.
    _objectify_bool = objectify.PyType('bool', None, objectify.BoolElement,
                                       str)
    _objectify_bool.register()
    # Serialize with 13 digits of precision.
    _objectify_float = objectify.PyType('float', None, objectify.FloatElement,
                                        lambda f: '%.13g' % f)
    _objectify_float.register()


def _unregister_objectify_types():
    """Clean up custom lxml.objectify types."""
    # The register/unregister process should really be performed as a context
    # manager to ensure the unregister happens even if an exception is raised in
    # between. A function annotation might be a good implementation option.
    global _objectify_bool, _objectify_float
    _objectify_bool.unregister()
    _objectify_float.unregister()


def _wellname_to_row(w):
    return ord(w[0]) - ord('A')


def _wellname_to_column(w):
    return int(w[1:]) - 1


def export_hpdd(design, agents, filename, assay_volume=60, dmso_limit=0.01,
                exclude_outer=True):
    """Create a .hpdd format file from a drug treatment design table.

    Parameters
    ----------
    design : DataFrame
        Must contain at least the following columns: well, agent, concentration.
        If the design contains more than one plate then a 'plate' column must
        also be present. The well column must contain A1-style well names. For
        drug combination designs, append more paired agent and concentration
        columns, e.g. 'agent2' and 'concentration2', then 'agent3' and
        'concentration3', and so on. Any suffix may be used as long as it is
        consistent between each agent/concentration column pair.
    agents : DataFrame
        Must contain at least the following columns: name, stock_concentration.
        The name column must match up with the values in the design.agent
        column (as well as any subsequent agent columns). All agents in the
        design table must be present in the agents table.
    filename : str
        Name of the file to write. Should end with '.hpdd'.
    assay_volume : Optional[float]
        Initial liquid volume present in each well before drug dispensing takes
        place, in microliters.
    dmso_limit: Optional[float]
        Percent of DMSO allowed in each well after drug dispensing.
    exclude_outer : bool
       exludes outer wells (wells that are not being used as negative control)
       from being back filled with DMSO.

    Notes
    -----
    All concentrations (design.concentration and agents.stock_concentration)
    must be in micromolar.

    Only 384 and 96 well plates are supported.

    The dmso_limit is not yet verified by this function. In the interim, the
    HPDD software may by used for this. (Open the generated file in the
    software and look for the yellow triangle warning symbol in any wells.)

    The hardware in the drug printer can only produce droplets of a specific
    volume. This puts a lower limit on the achievable drug concentrations as
    well as constrains the specific concentrations achievable at the lower end
    of the range. This function does not implement checks on these limitations
    yet, but again the HPDD software may be used for verification.

    """
    if 'plate' not in design:
        if 'barcode' in design:
            design['plate'] = design['barcode']
        else:
            design = design.copy()
            design['plate'] = 'plate_1'
    plate_names = design.plate.unique()

    plate_size = len(design) / len(plate_names)
    if plate_size not in (96, 384):
        raise ValueError("Only complete 96 and 384 well plates are supported")
    # Math to calculate number of rows and columns from total well count. The
    # key value (called `plate_scale` here) is the base 4 log of the plate size
    # divided by 6 -- log_4(s/6).
    plate_scale = np.log(plate_size / 6) / np.log(4)
    plate_columns = int(3 * (2 ** plate_scale))
    plate_rows = int(2 * (2 ** plate_scale))

    agent_columns = design.columns[design.columns.str.startswith('agent')]
    agent_column_suffixes = agent_columns.str.replace('agent', '')
    design_agents = [a for a in design[agent_columns].stack().unique()
                     if a not in  ['', 'DMSO']]
    if sorted(design_agents) != sorted(agents.name.unique()):
        raise ValueError("Mismatch between design and agents tables")

    _register_objectify_types()

    root = objectify.E.Protocol()
    created_comment = 'created by DataRail version %s' % package_version
    root.append(etree.Comment(' %s ' % created_comment))
    root.Version = 3
    root.VolumeUnit = 'nL'
    root.ConcentrationUnit = MU + 'M'
    root.MolarityConcentrationUnit = MU + 'M';
    root.MassConcentrationUnit = 'ng_mL'

    root.Fluids = objectify.E.Fluids()

    root.ShakePerFluid = False
    root.ShakePlateDuration = 5
    root.ShakePerWell = True
    root.ShakeThresholdVolume = 100
    root.BackfillOrder = 'Last'
    root.BackfillNoDispense = False

    root.Plates = objectify.E.Plates()
    root.Backfills = objectify.E.Backfills()

    fluid_ids = {}
    for a in agents.itertuples():
        fluid = objectify.E.Fluid()
        fluid_id = str(a.Index)
        fluid.set('ID', fluid_id)
        fluid.Name = str(a.name)
        fluid.ByVol = False
        # This value is relative to the top-level ConcentrationUnit.
        fluid.Concentration = a.stock_concentration
        # This is the unit for display in the HPDD software only!
        fluid.ConcentrationUnit = 'mM'
        fluid.ClassID = 0
        root.Fluids.append(fluid)
        fluid_ids[a.name] = fluid_id

    plate_ids = {}
    for pi, p in enumerate(plate_names):
        plate_ids[p] = pi
        plate = objectify.E.Plate()
        plate.PlateType = 'Default%d' % plate_size
        plate.Rows = plate_rows
        plate.Cols = plate_columns
        plate.Name = str(p)
        plate.AssayVolume = int(assay_volume * 1e3)
        plate.DMSOLimit = dmso_limit
        plate.DontShake = False
        plate.Wells = objectify.E.Wells()
        root.Plates.append(plate)
        for d in design.loc[design.plate == p].itertuples():
            well = objectify.E.Well()
            well.set('Row', str(_wellname_to_row(d.well)))
            well.set('Col', str(_wellname_to_column(d.well)))
            for a_column, suffix in zip(agent_columns, agent_column_suffixes):
                c_column = 'concentration' + suffix
                # TODO: Round concentration based on D300 droplet volumes and
                # emit a warning if there is a significant discrepancy.
                # TODO: Warn if requested DMSO limit is violated.
                agent_concentration = float(getattr(d, c_column))
                agent_name = getattr(d, a_column)
                if ((np.isreal(agent_name) and np.isnan(agent_name))
                    or agent_name == ''):
                    continue
                well_fluid = objectify.E.Fluid(agent_concentration)
                well_fluid.set('ID', fluid_ids[agent_name])
                well.append(well_fluid)
            if well.countchildren():
                plate.Wells.append(well)
        plate.Randomize = objectify.E.Randomize()

    backfill = objectify.E.Backfill()
    backfill.set('ClassID', str(0))
    backfill.set('Type', 'ToMaxVolume')
    root.Backfills.append(backfill)
    backfill_wells = objectify.E.Wells()
    backfill.append(backfill_wells)
    if exclude_outer:
        design = design.copy()
        design = design[~design['role'].isnull()]
    for p in plate_names:
        pi = plate_ids[p]
        for d in design.loc[design.plate == p].itertuples():
            well = objectify.E.Well()
            well.set('P', str(pi))
            well.set('R', str(_wellname_to_row(d.well)))
            well.set('C', str(_wellname_to_column(d.well)))
            backfill_wells.append(well)

    objectify.deannotate(root, cleanup_namespaces=True)
    content = etree.tostring(root, pretty_print=True, xml_declaration=True,
                             encoding='utf-8')

    _unregister_objectify_types()

    with open(filename, 'wb') as f:
        f.write(content)


def import_hpdd(filename):
    raise NotImplementedError()
