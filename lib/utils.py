from pathlib import Path


def write_manual_pe_scheme(
    outfile: str | Path,
    pe_dir: str,
    n_volumes: int,
    readout_time: float = 0.1,
):
    """
    Reproduces dwifslpreproc behavior when:
    - pe_dir is provided
    - readout_time is NOT provided (defaults to 0.1)
    """
    from mrtrix3 import phaseencoding

    pe_vec = list(phaseencoding.direction(pe_dir))    # [dx, dy, dz]

    with open(outfile, "w") as f:
        for _ in range(n_volumes):
            f.write(f"{pe_vec[0]} {pe_vec[1]} {pe_vec[2]} {readout_time}\n")
