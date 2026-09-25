"""Parameter sections, presets and project.yaml round trips."""

import pytest

from atlas_ms.config import AdductSettings, InstrumentSettings, ProjectConfig


def test_yaml_round_trip(tmp_path):
    config = ProjectConfig()
    config.instrument.noise_threshold = 5e3
    config.save(tmp_path / "project.yaml")
    loaded = ProjectConfig.load(tmp_path / "project.yaml")
    assert loaded.to_dict() == config.to_dict()


def test_missing_values_take_defaults():
    config = ProjectConfig.from_dict({"linking": {"rt_tol_s": 12.0}})
    assert config.linking.rt_tol_s == 12.0
    assert config.linking.mz_tol_ppm == ProjectConfig().linking.mz_tol_ppm


def test_unknown_keys_are_rejected():
    with pytest.raises(ValueError, match="Unknown parameter"):
        ProjectConfig.from_dict({"instrument": {"noise_treshold": 1e4}})  # typo
    with pytest.raises(ValueError, match="Unknown section"):
        ProjectConfig.from_dict({"instrumnet": {}})


def test_values_are_validated():
    with pytest.raises(ValueError):
        InstrumentSettings(noise_threshold="high")
    with pytest.raises(ValueError):
        InstrumentSettings(mass_error_ppm=-1.0)


def test_instrument_preset_sets_several_sections():
    config = ProjectConfig()
    config.apply_instrument_preset("qtof")
    assert config.presets.instrument == "qtof"
    assert config.instrument.mass_error_ppm == 20.0
    assert config.linking.mz_tol_ppm == 15.0
    with pytest.raises(ValueError, match="Unknown instrument preset"):
        config.apply_instrument_preset("fticr")


def test_openms_adducts_are_normalised():
    adducts = AdductSettings(adducts=[
        {"name": "[M+H]+", "openms": "H:+", "probability": 2.0},
        {"name": "[M+Na]+", "openms": "Na:+", "probability": 2.0},
        {"name": "-H2O", "openms": "H-2O-1:0", "probability": 0.05},
    ])
    assert adducts.openms_adducts() == [b"H:+:0.5", b"Na:+:0.5", b"H-2O-1:0:0.05"]


def test_adduct_polarity_is_checked():
    with pytest.raises(ValueError, match="does not match positive"):
        AdductSettings(polarity="positive", adducts=[{"name": "[M-H]-", "openms": "H-1:-", "probability": 1}])
    config = ProjectConfig()
    config.apply_adduct_preset("negative_default")
    assert config.adducts.polarity == "negative"
