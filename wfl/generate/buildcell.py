import io
import numbers
import subprocess

import ase.data
import numpy as np
import spglib
from ase.atoms import Atoms
from ase.io.extxyz import key_val_str_to_dict

from wfl.autoparallelize import _autoparallelize_ll
from wfl.utils.round_sig_figs import round_sig_figs


# might be nice to combine Z, vol_per_atom, and bond_lengths into a single dict with Z as key and vol
# and bond len as values
def create_input(z, vol_per_atom, bond_lengths, filename, composition=None, RSS_min_vol_factor=0.5,
                 vol_range=(0.95, 1.05), min_sep_factor=0.9, symmops='1-8', natom=(6, 24), odd=None, verbose=False):
    if composition is None:
        composition = [1]
    if isinstance(z, numbers.Integral):
        z = [z]
    if isinstance(vol_per_atom, numbers.Real):
        vol_per_atom = [vol_per_atom]
    if isinstance(bond_lengths, numbers.Real):
        bond_lengths = [bond_lengths]

    if len(z) != len(composition) or len(z) != len(vol_per_atom) or len(z) != len(bond_lengths):
        raise ValueError(
            'Got mismatched lengths for Z {} composition {} vol_per_atom {} minsep {}'.format(len(z), len(composition),
                                                                                              len(vol_per_atom),
                                                                                              len(bond_lengths)))

    assert odd is None or odd == 'only' or odd == 'also'

    min_sep_str = '0.5'
    for (Z_1, bond_len_1) in zip(z, bond_lengths):
        for (Z_2, bond_len_2) in zip(z, bond_lengths):
            if Z_1 <= Z_2:
                min_sep_str += ' {}-{}={}'.format(ase.data.chemical_symbols[Z_1], ase.data.chemical_symbols[Z_2],
                                                  round_sig_figs(min_sep_factor * 0.5 * (bond_len_1 + bond_len_2), 2))
    if verbose:
        print('Z', z, 'vol_per_atom', vol_per_atom, 'composition', composition)
    species_str = ','.join(['{}%NUM={}'.format(ase.data.chemical_symbols[z], n) for (z, n) in zip(z, composition)])
    vol_per_formula_unit = np.sum([v * n for (v, n) in zip(vol_per_atom, composition)])
    vol_per_atom = vol_per_formula_unit / sum(composition)
    target_volume = vol_per_atom * np.count_nonzero(np.array(composition) != 0)

    nform = []
    nat_per_fu = np.sum(composition)
    for nat in range(natom[0], natom[1] + 1):
        if nat % nat_per_fu != 0:  # incompatible with composition
            continue
        if nat % 2 == 0 and odd == 'only':  # even but odd_only was specified
            continue
        if nat % 2 == 1 and odd is None:  # odd
            continue
        nform.append(int(nat / nat_per_fu))
    nform_str = '{' + ','.join([str(n) for n in nform]) + '}'

    with open(filename, 'w') as fout:
        # fixme: add docstring for the craziness of buildcell here
        fout.write('#TARGVOL={}-{}\n'.format(round_sig_figs(target_volume * vol_range[0], 2),
                                             round_sig_figs(target_volume * vol_range[1], 2)))
        fout.write('#SPECIES={}\n'.format(species_str))
        fout.write('#NFORM={}\n'.format(nform_str))
        fout.write('#SYMMOPS={}\n'.format(symmops))
        fout.write('#SLACK=0.25\n')
        fout.write('#OVERLAP=0.1\n')
        fout.write('#COMPACT\n')
        fout.write('#MINSEP={}\n'.format(min_sep_str))
        fout.write('##EXTRA_INFO RSS_min_vol_per_atom={}\n'.format(RSS_min_vol_factor * vol_per_atom))


# could this be replaced with function from ase.io.castep without too much loss of speed?
def conv_buildcell_out(buildcell_output):
    """Convert from buildcell output to atoms object(s)

    Parameters
    ----------
    buildcell_output: str
        stdout of buildcell run

    Returns
    -------
    ats: list(Atoms)
        list of created Atoms objects
    """
    output_file = io.StringIO(buildcell_output)

    cell = None
    species = None
    pos_frac = None

    ats = []
    in_cell = False
    in_pos = False
    for l in output_file:
        if in_cell:
            cell = [float(a) for a in l.split()]
            l = next(output_file)
            cell.extend([float(a) for a in l.split()])
            in_cell = False
            pos_frac = []
            species = []
        if in_pos:
            if '%ENDBLOCK POSITIONS_FRAC' in l:
                in_pos = False
                if cell is None or species is None or pos_frac is None:
                    raise RuntimeError(('got to ENDBLOCK POSITIONS_FRAC without one of cell, species, '
                                        'pos_fract (found? {} {} {}, respectively)').format(cell is None,
                                                                                            species is None,
                                                                                            pos_frac is None))
                at = Atoms(cell=cell, symbols=species, scaled_positions=pos_frac, pbc=[True] * 3)
                ats.append(at)
                # ase.io.write(sys.stdout, at, format='extxyz')
            else:
                f = l.split()
                species.append(f[0])
                pos_frac.append([float(x) for x in f[1:]])
        if '%BLOCK LATTICE_ABC' in l:
            in_cell = True
        if '%BLOCK POSITIONS_FRAC' in l:
            in_pos = True

    return ats


def run(outputs, config_is, buildcell_cmd, buildcell_input, extra_info=None,
        perturbation=0.0, skip_failures=True, symprec=0.01, verbose=False):
    """Creates atomic configurations by repeatedly running buildcell, I/O with OutputSpec

    Parameters
    ----------
    outputs: OutputSpec
        where to write outputs
    config_is: int / list(int)
        numbers to set in buildcell_config_i info field
        length sets number of configs generated
    buildcell_cmd:
        `buildcell` executable, including path if needed
    buildcell_input: str
        input file contents to buildcell, in the form of a single string
    extra_info: dict, default {}
        extra fields to place into created atoms info dict
    perturbation: float, default 0.0
        magnitude of perturbation to atomic positions with Atoms.rattle()
    skip_failures: bool, default True
        skip failed buildcell calls
    symprec: float, default 0.01
        precision for symmetry check
    verbose: bool, default False
        print some running info

    Returns
    ------
        ConfigSet corresponding to output
    """
    if extra_info is None:
        extra_info = {}
    return _autoparallelize_ll(iterable=config_is, outputspec=outputs, op=run_autopara_wrappable,
                         buildcell_cmd=buildcell_cmd, buildcell_input=buildcell_input,
                         extra_info=extra_info, perturbation=perturbation, skip_failures=skip_failures,
                         verbose=verbose, symprec=symprec)


def run_autopara_wrappable(config_is, buildcell_cmd, buildcell_input, extra_info=None,
           perturbation=0.0, skip_failures=True, symprec=0.01, verbose=False):
    """Creates atomic configurations by repeatedly running buildcell

    Parameters
    ----------
    config_is: int / list(int)
        numbers to set in buildcell_config_i info field
        length sets number of configs generated
    buildcell_cmd:
        `buildcell` executable, including path if needed
    buildcell_input: str
        input file contents to buildcell, in the form of a single string
    extra_info: dict, default {}
        extra fields to place into created atoms info dict
    perturbation: float, default 0.0
        magnitude of perturbation to atomic positions with Atoms.rattle()
    skip_failures: bool, default True
        skip failed buildcell calls
    symprec: float, default 0.01
        precision for symmetry check
    verbose: bool, default False
        print some running info

    Returns
    ------
        list of Atoms created by buildcell
    """

    if extra_info is None:
        extra_info = {}
    if isinstance(config_is, int):
        config_is = [config_is]

    if verbose:
        print('repeat_buildcell calling with stdin ' + buildcell_input)

    t_extra_info = {}
    for l in buildcell_input.splitlines():
        if '##EXTRA_INFO' in l:
            t_extra_info.update(key_val_str_to_dict(l.replace('##EXTRA_INFO', '').strip()))
    t_extra_info.update(extra_info)
    extra_info = t_extra_info

    atoms_list = []
    for config_i in config_is:
        p = subprocess.run(buildcell_cmd, input=buildcell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           encoding='utf-8')
        if len(p.stdout) == 0:
            raise RuntimeError('Got no output from buildcell')
        if verbose:
            print('repeat_buildcell got stdout' + p.stdout)

        ats = None
        try:
            ats = conv_buildcell_out(p.stdout)
        except Exception as exc:
            if skip_failures:
                print(f'buildcell returned bad output exception {exc}')
            else:
                raise exc

        if ats is not None:
            if len(ats) != 1:
                raise RuntimeError(
                    'Got unexpected number of atoms {} != 1 from converting buildcell output'.format(len(ats)))
            at0 = ats[0]
            at0.info['config_type'] = 'buildcell'
            at0.info['buildcell_config_i'] = config_i
            at0.info.update(extra_info)
            at0.rattle(perturbation)

            # handle the symmetry
            dataset = spglib.get_symmetry_dataset((at0.cell, at0.get_scaled_positions(), at0.numbers), symprec=symprec)
            at0.info['buildcell_symmetry'] = '{} {} {} @ {}'.format(dataset['number'], dataset['international'],
                                                                    dataset['hall'], symprec)
            atoms_list.append(at0)

    return atoms_list
