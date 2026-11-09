"""Normalize scanner-exported gradient coefficient files for the GNC backends."""

from pathlib import Path


def strip_coeff_header(coeff_file: str | Path, scratch_dir: str | Path) -> str:
    """
    Strip '#'-prefixed header comments from a scanner-exported gradient
    coefficient file, producing the format gradient_unwarp.py and
    TORTOISE's CreateGradientNonlinearityBMatrix expect.

    Args:
        coeff_file: Path to the scanner gradient coefficient file
        scratch_dir: Directory to write the stripped copy into, if needed

    Returns:
        coeff_file unchanged if it has no '#' lines, otherwise the path to
        a stripped copy written into scratch_dir
    """
    # newline='' disables universal-newline translation so CRLF line endings
    # (as in real scanner exports) round-trip unchanged.
    with open(coeff_file, newline='') as f:
        lines = f.readlines()
    stripped = [line for line in lines if not line.lstrip().startswith('#')]
    if len(stripped) == len(lines):
        return str(coeff_file)

    out_path = Path(scratch_dir) / f'{Path(coeff_file).stem}_stripped.grad'
    with open(out_path, 'w', newline='') as f:
        f.writelines(stripped)
    return str(out_path)
