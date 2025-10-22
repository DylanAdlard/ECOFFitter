import pytest
import numpy as np
import pandas as pd
from scipy.optimize import OptimizeResult
from ecoff_fitter import EcoffFitter


@pytest.fixture
def wt_samples():
    """Fixture for sample data (wild-type, no mutations)."""
    samples = pd.DataFrame({"UNIQUEID": ["A", "B", "C"], "MIC": ["1", "2", "3"]})
    mutations = pd.DataFrame({"UNIQUEID": ["A"], "MUTATION": [None]})
    return samples, mutations


@pytest.fixture
def mixed_variants():
    """Fixture for data with both wild-type and mutant samples."""
    samples = pd.DataFrame(
        {"UNIQUEID": ["A", "B", "C", "D"], "MIC": ["1", "<=2", ">3", "4"]}
    )
    mutations = pd.DataFrame(
        {
            "UNIQUEID": ["A", "B", "C"],
            # B has a synonymous change, C has a non-synonymous change
            "MUTATION": [None, "mut1@G12G", "mut2@A13V"],
        }
    )
    return samples, mutations


def test_flag_test1_wt(wt_samples, mixed_variants):
    # Test gWT_definition="test1" filtering
    samples, mutations = wt_samples
    ecoff = EcoffFitter(samples, mutations, gWT_definition="test1")
    # All samples should be flagged as WT
    assert all(ecoff.df["WT"]), f"Expected all WT flags True, got {list(ecoff.df['WT'])}"

    samples, mutations = mixed_variants
    ecoff = EcoffFitter(samples, mutations, gWT_definition="test1")
    # A: no mutation -> WT; B: synonymous -> WT; C: non-synonymous -> not WT; D: no entry -> WT
    expected_wt_flags = [True, True, False, True]
    assert list(ecoff.df["WT"]) == expected_wt_flags, (
        f"Expected WT flags {expected_wt_flags}, got {list(ecoff.df['WT'])}"
    )


def test_define_intervals_uncensored(wt_samples):
    # Uncensored interval definition uses tail_dilutions
    samples, mutations = wt_samples
    ecoff = EcoffFitter(samples, mutations, censored=False, tail_dilutions=1)
    y_low, y_high = ecoff.define_intervals()
    # Exact values: [1/2, 2/2, 3/2] and [1,2,3]
    expected_low = np.array([0.5, 1.0, 1.5])
    expected_high = np.array([1.0, 2.0, 3.0])
    # Log2 transform
    log2 = lambda x: np.log(x) / np.log(2)
    assert np.allclose(y_low, log2(expected_low), atol=1e-3)
    assert np.allclose(y_high, log2(expected_high), atol=1e-3)


def test_define_intervals_censored(mixed_variants):
    # Censored interval definition for mixed variants
    samples, mutations = mixed_variants
    ecoff = EcoffFitter(samples, mutations, censored=True)
    y_low, y_high = ecoff.define_intervals()
    # Expected: 
    # A: [0.5,1] => log2: [-1,0]
    # B: left-censored <=2 => [1e-6,2] => [log2(1e-6),1]
    # C: right-censored >3 => [3,inf] => [log2(3), inf]
    # D: [2,4] => [1,2]
    assert pytest.approx(-1.0, abs=1e-3) == y_low[0]
    assert pytest.approx(0.0, abs=1e-3) == y_high[0]
    assert pytest.approx(np.log(1e-6)/np.log(2), abs=1e-3) == y_low[1]
    assert pytest.approx(1.0, abs=1e-3) == y_high[1]
    assert pytest.approx(np.log(3)/np.log(2), abs=1e-3) == y_low[2]
    assert y_high[2] == np.inf
    assert pytest.approx(1.0, abs=1e-3) == y_low[3]
    assert pytest.approx(2.0, abs=1e-3) == y_high[3]


def test_log_transf_intervals(wt_samples):
    # Test direct log transformation
    samples, mutations = wt_samples
    ecoff = EcoffFitter(samples, mutations)
    y_low = np.array([0.5, 1.0, 1.5])
    y_high = np.array([1.0, 2.0, 3.0])
    y_low_log, y_high_log = ecoff.log_transf_intervals(y_low, y_high)
    expected_low = np.array([np.log(0.5)/np.log(2), 0.0, np.log(1.5)/np.log(2)])
    expected_high = np.array([0.0, 1.0, np.log(3)/np.log(2)])
    assert np.allclose(y_low_log, expected_low, atol=1e-3)
    assert np.allclose(y_high_log, expected_high, atol=1e-3)


def test_fit_model_returns_optimize_result(wt_samples):
    samples, mutations = wt_samples
    ecoff = EcoffFitter(samples, mutations)
    result = ecoff.fit()
    assert isinstance(result, OptimizeResult)
    assert np.isfinite(result.x[0])
    assert np.isfinite(result.x[1])


def test_generate_ecoff_basic(wt_samples):
    samples, mutations = wt_samples
    ecoff = EcoffFitter(samples, mutations)
    ecoff_value, z_percentile, mu, sigma, model = ecoff.generate(percentile=99)
    # Basic sanity checks
    assert isinstance(ecoff_value, float)
    assert isinstance(model, OptimizeResult)
    assert z_percentile > mu
    assert sigma > 0
    # ECOFF should equal dilution_factor**z_percentile
    assert ecoff_value == pytest.approx(
        ecoff.dilution_factor ** z_percentile, rel=1e-6
    )
