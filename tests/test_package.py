"""Basic package tests to verify installation."""

import kappenball


def test_package_import():
    """Test that kappenball package can be imported."""
    assert kappenball.__version__ == "0.1.0"


def test_modules_exist():
    """Test that all expected modules are available."""
    assert hasattr(kappenball, "billiards")
    assert hasattr(kappenball, "falling_ball")
    assert hasattr(kappenball, "physics")
    assert hasattr(kappenball, "visualisation")


def test_module_imports():
    """Test that modules can be imported directly."""
    from kappenball import billiards, falling_ball, physics, visualisation
    
    assert billiards is not None
    assert falling_ball is not None
    assert physics is not None
    assert visualisation is not None

