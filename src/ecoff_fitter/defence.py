import pandas as pd

def validate_input(
    samples, mutations, gWT_definition, dilution_factor, censored, tail_dilutions
):
    """Validates inputs for the ECOFF generator initialization."""

    assert isinstance(samples, pd.DataFrame), "samples must be a pandas DataFrame."

    # Check required columns in samples
    assert all(
        column in samples.columns for column in ["UNIQUEID", "MIC"]
    ), "Input samples must contain columns 'UNIQUEID' and 'MIC'"

    if gWT_definition is not None:
        assert isinstance(mutations, pd.DataFrame), "mutations must be a pandas DataFrame."
        assert gWT_definition in ['test1', 'ERJ2022'], 'only test1 and ERJ2022 gWT protocols are implemented'
        assert all(
            column in mutations.columns for column in ["UNIQUEID", "MUTATION"]
        ), "Input mutations must contain columns 'UNIQUEID' and 'MUTATION'"

    # Validate dilution_factor
    assert (
        isinstance(dilution_factor, int) and dilution_factor > 0
    ), "dilution_factor must be a positive integer."

    # Validate censored flag
    assert isinstance(
        censored, bool
    ), "censored must be a boolean value (True or False)."

    # Validate tail_dilutions if censored is False
    if not censored:
        assert (
            isinstance(tail_dilutions, int) and tail_dilutions > 0
        ), "When censored is False, tail_dilutions must be a positive integer or specified."



