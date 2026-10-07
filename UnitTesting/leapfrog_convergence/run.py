#!/usr/bin/env python3
"""Run real C convergence regressions with Python 3 and a C99 compiler.

    python3 UnitTesting/leapfrog_convergence/run.py
    CC=clang CFLAGS='-O1 -g -fsanitize=address,undefined' python3 UnitTesting/leapfrog_convergence/run.py

CC (default: cc), CFLAGS (default: -O2), and LDFLAGS are respected. Builds
and executables live in a temporary directory. --source-root can target an
unpatched checkout to verify that the regressions fail there. These tests
stub only the Cactus interface; full Einstein Toolkit testing is still needed.
"""

import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


ROUTINES = """
MoL_malloc_y_n_gfs MoL_malloc_non_y_n_gfs MoL_free_memory_y_n_gfs
MoL_free_memory_non_y_n_gfs MoL_step_forward_in_time set_Cparameters_to_default
set_Nxx_dxx_invdx_params__and__xx driver_bcstruct set_up__bc_gz_map_and_parity_condns
set_bcstruct freemem_bcstruct find_timestep rfm_precompute_rfmstruct_malloc
rfm_precompute_rfmstruct_define rfm_precompute_rfmstruct_freemem
initial_guess_all_points initial_guess_single_point auxevol_gfs_all_points
wavespeed_gf_all_points rhs_eval residual_all_points L2_norm_of_gf L2_norm_weights
leapfrog_step_forward_in_time set_optimum_eta_damping
""".split()


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=here.parents[1])
    args = parser.parse_args()
    src = args.source_root.resolve() / "NRPyEllipticET/src/conformally_flat_BBH"
    files = [src / ("conformally_flat_BBH_" + name + ".c") for name in ROUTINES]
    files += [src / ("NRPyEllipticET_conformally_flat_BBH_" + name + ".c")
              for name in ("sequencing", "Hyperbolic_Relaxation")]
    failures = 0
    with tempfile.TemporaryDirectory(prefix="nrpy-convergence-") as tmp:
        executable = str(Path(tmp) / "check")
        command = shlex.split(os.environ.get("CC", "cc"))
        command += ["-std=c99", *shlex.split(os.environ.get("CFLAGS", "-O2")),
                    "-I" + str(here), "-I" + str(src), "-I" + str(src.parent)]
        command += [str(path) for path in files] + [str(here / "check.c")]
        command += [*shlex.split(os.environ.get("LDFLAGS", "")), "-lm", "-o", executable]
        subprocess.run(command, check=True)
        for case in ("sequencing-cached", "standard-cached", "sequencing-final", "standard-exhausted"):
            result = subprocess.run([executable, case], text=True, capture_output=True)
            expected = 2 if case == "standard-exhausted" else 0
            passed = result.returncode == expected
            if case == "standard-exhausted":
                passed = passed and "Solution did not reach the specified tolerance" in result.stderr
            print(("PASS " if passed else "FAIL ") + case, flush=True)
            if result.stdout:
                print(result.stdout, end="", flush=True)
            if not passed:
                print("exit status:", result.returncode)
                print(result.stderr, end="", flush=True)
                failures += 1
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
