"""
VASP calculator

"""

import os

import json

import numpy as np

from ase.calculators.calculator import all_changes
from ase.calculators.vasp.vasp import Vasp as ASE_Vasp

from .wfl_fileio_calculator import WFLFileIOCalculator
from .utils import handle_nonperiodic

# NOMAD compatible, see https://nomad-lab.eu/prod/rae/gui/uploads
_default_keep_files = ["POSCAR", "INCAR", "KPOINTS", "OUTCAR", "vasprun.xml", "vasp.out"]
_default_properties = ["energy", "forces", "stress"]

class Vasp(WFLFileIOCalculator, ASE_Vasp):
    """Extension of ASE's Vasp calculator that can be used by wfl.calculators.generic

    Notes
    -----
    "directory" argument cannot be present. Use rundir_prefix and workdir instead.
    "pp" defaults to ".", so VASP_PP_PATH env var is absolute path to "<elem name>/POTCAR" files

    Parameters
    ----------
    keep_files: bool / None / "default" / list(str), default "default"
        what kind of files to keep from the run
            True, "*" : everything kept
            None, False : nothing kept
            "default"   : only ones needed for NOMAD uploads ('POSCAR', 'INCAR', 'KPOINTS', 'OUTCAR', 'vasprun.xml', 'vasp.out')
            list(str)   : list of file globs to save
    rundir_prefix: str / Path, default 'run\_VASP\_'
        Run directory name prefix
    workdir: str / Path, default . at calculate time
        Path in which rundir will be created.
    scratchdir: str / Path, default None
        temporary directory to execute calculations in and delete or copy back results (set by
        "keep_files") if needed.  For example, directory on a local disk with fast file I/O.
    calculator_exec: str
        executable to run (without ASE-specific command line arguments). Mutually exclusive with ASE-built-in "command"
    calculator_exec_gamma: str
        executable to run for nonperiodic systems (overrides ASE_VASP_COMMAND_GAMMA, VASP_COMMAND_GAMMA, VASP_SCRIPT_GAMMA).
        Mutually exclusive with ASE-built-in "command" and "command_gamma"
    **kwargs: arguments for ase.calculators.vasp.vasp.Vasp
        remaining arguments to ASE's Vasp calculator constructor
    """

    # default value of wfl_num_inputs_per_python_subprocess for calculators.generic,
    # to override that function's built-in default of 10
    wfl_generic_num_inputs_per_python_subprocess = 1

    # note that we also have a default for "pp", but that has to be handlded separately
    default_parameters = ASE_Vasp.default_parameters.copy()
    default_parameters.update({
        "ismear": 0,
        "lwave": False,
        "lcharg": False,
    })

    def __init__(self, keep_files="default", rundir_prefix="run_VASP_",
                 workdir=None, scratchdir=None,
                 calculator_exec=None, calculator_exec_gamma=None,
                 **kwargs):

        # get initialparams from env var
        kwargs_use = {}
        if 'WFL_VASP_KWARGS' in os.environ:
            try:
                kwargs_use = json.loads(os.environ['WFL_VASP_KWARGS'])
            except json.decoder.JSONDecodeError:
                with open(os.environ['WFL_VASP_KWARGS']) as fin:
                    kwargs_use = json.load(fin)

        # override with explicitly passed in values
        kwargs_use.update(kwargs)

        # pp is not handled by Vasp.default_parameters, because that only includes
        # parameters that are in INCAR
        if "pp" not in kwargs_use:
            kwargs_use["pp"] = "."

        if calculator_exec is not None:
            if "command" in kwargs_use or "command_gamma" in kwargs_use:
                raise ValueError("Cannot specify both calculator_exec and command or command_gamma")
            self.command = calculator_exec
        if calculator_exec_gamma is not None:
            if "command" in kwargs_use or "command_gamma" in kwargs_use:
                raise ValueError("Cannot specify both calculator_exec_gamma and command or command_gamma")
            self._command_gamma = calculator_exec_gamma
        else:
            self._command_gamma = kwargs_use.pop("command_gamma", None)

        # WFLFileIOCalculator is a mixin, will call remaining superclass constructors for us
        super().__init__(keep_files=keep_files, rundir_prefix=rundir_prefix,
                         workdir=workdir, scratchdir=scratchdir, **kwargs_use)

    def per_config_setup(self, atoms, nonperiodic):
        # VASP cannot handle nonperiodic, and current Vasp calculator complains if pbc=F
        self._orig_pbc = atoms.pbc.copy()
        atoms.pbc = [True] * 3

        # VASP requires periodic cells with non-zero cell vectors
        if np.dot(np.cross(atoms.cell[0], atoms.cell[1]), atoms.cell[2]) < 0.0:
            t = atoms.cell[0].copy()
            atoms.cell[0] = atoms.cell[1]
            atoms.cell[1] = t
            self._permuted_a0_a1 = True
        else:
            self._permuted_a0_a1 = False

        self._orig_command = self.command
        if nonperiodic:
            if self._command_gamma is not None:
                command_gamma = self._command_gamma
            else:
                command_gamma = None
                for env_var in self.env_commands:
                    if env_var + "_GAMMA" in os.environ:
                        command_gamma = os.environ[env_var + "_GAMMA"]
                        break
            if command_gamma is not None:
                self.command = command_gamma

        # set some things for nonperiodic systems
        self._orig_kspacing = self.float_params["kspacing"]
        self._orig_kgamma = self.bool_params["kgamma"]
        if nonperiodic:
            self.float_params["kspacing"] = 100000.0
            self.bool_params["kgamma"] = True


    def per_config_restore(self, atoms):
        # undo pbc change
        atoms.pbc[:] = self._orig_pbc

        # undo cell vector permutation
        if self._permuted_a0_a1:
            t = atoms.cell[0].copy()
            atoms.cell[0] = atoms.cell[1]
            atoms.cell[1] = t

        # undo nonperiodic related changes
        self.command = self._orig_command
        self.float_params["kspacing"] = self._orig_kspacing
        self.bool_params["kgamma"] = self._orig_kgamma


    def calculate(self, atoms=None, properties=_default_properties, system_changes=all_changes):
        """Do the calculation. Handles the working directories in addition to regular
        ASE calculation operations (writing input, executing, reading_results)"""

        if atoms is not None:
            self.atoms = atoms.copy()

        # hack to disable parallel I/O, since that makes VASP's I/O break if
        # this script is being run with mpirun
        from ase.parallel import world, DummyMPI
        orig_world_comm = world.comm
        world.comm = DummyMPI()

        # make sure VASP_PP_PATH is set to a sensible default
        orig_VASP_PP_PATH = os.environ.get("VASP_PP_PATH")
        if orig_VASP_PP_PATH is None:
            if self.input_params['pp'].startswith("/"):
                os.environ["VASP_PP_PATH"] = "/."
            else:
                os.environ["VASP_PP_PATH"] = "."

        nonperiodic, properties_use = handle_nonperiodic(atoms, properties)
        # from WFLFileIOCalculator
        self.per_config_setup(atoms, nonperiodic)

        # from WFLFileIOCalculator
        self.setup_rundir()

        atoms.info["vasp_rundir"] = str(self._cur_rundir)

        if self.float_params["encut"] is None:
            raise RuntimeError("Refusing to run without explicit ENCUT")
        # should we require anything except ENCUT?

        calculation_succeeded = False
        try:
            super().calculate(atoms=atoms, properties=properties_use, system_changes=system_changes)
            calculation_succeeded = True
            if 'DFT_FAILED_VASP' in atoms.info:
                del atoms.info['DFT_FAILED_VASP']
        except Exception as exc:
            atoms.info['DFT_FAILED_VASP'] = True
            raise exc
        finally:
            # from WFLFileIOCalculator
            self.clean_rundir(_default_keep_files, calculation_succeeded)

            self.per_config_restore(atoms)

            # undo env VASP_PP_PATH
            if orig_VASP_PP_PATH is not None:
                os.environ["VASP_PP_PATH"] = orig_VASP_PP_PATH

            # undo communicator mangling
            world.comm = orig_world_comm
