"""Tests for pytheas -- physics validation and regression."""

import numpy as np
import pytest
from dataclasses import asdict, FrozenInstanceError
from datetime import datetime, timedelta

from pytheas import (
    compute_g,
    compute_timeseries,
    GravityResult,
    TimeSeries,
    normal_gravity,
    sun_position_ecef,
    moon_position_ecef,
    geodetic_to_ecef,
    enu_basis,
    measurement_axis,
    tidal_acceleration,
    julian_date,
    gmst_rad,
    GravityField,
    LabFrame,
    GM_MOON,
    GM_SUN,
    A_WGS84,
    OMEGA,
    DELTA_GRAV,
    GAMMA_E,
    GAMMA_P,
)


# =========================================================================
# Constants
# =========================================================================

class TestConstants:
    def test_gm_moon(self):
        assert abs(GM_MOON - 4.9028695e12) / 4.9028695e12 < 1e-6

    def test_gm_sun(self):
        assert abs(GM_SUN - 1.32712440018e20) / 1.32712440018e20 < 1e-6

    def test_gravimetric_factor(self):
        assert abs(DELTA_GRAV - 1.1608) < 0.001

    def test_equatorial_gravity(self):
        assert abs(GAMMA_E - 9.7803253141) < 1e-8

    def test_polar_gravity(self):
        assert abs(GAMMA_P - 9.8321849378) < 1e-8


# =========================================================================
# Normal gravity
# =========================================================================

class TestNormalGravity:
    def test_equator_sea_level(self):
        g = normal_gravity(0.0, 0.0)
        assert abs(g - 9.7803) < 0.001

    def test_pole_sea_level(self):
        g = normal_gravity(90.0, 0.0)
        assert abs(g - 9.8322) < 0.001

    def test_latitude_monotonic(self):
        """g increases with latitude (equator to pole)."""
        lats = [0, 15, 30, 45, 60, 75, 90]
        gs = [normal_gravity(lat, 0.0) for lat in lats]
        for i in range(len(gs) - 1):
            assert gs[i + 1] > gs[i]

    def test_altitude_decreases(self):
        """g decreases with altitude."""
        g0 = normal_gravity(45.0, 0.0)
        g1 = normal_gravity(45.0, 1000.0)
        assert g1 < g0

    def test_free_air_gradient(self):
        """Free-air gradient ~-3.086 uGal/m near sea level."""
        g0 = normal_gravity(45.0, 0.0)
        g1 = normal_gravity(45.0, 1.0)
        gradient = (g1 - g0) * 1e6  # um/s^2 per meter
        assert abs(gradient - (-3.086)) < 0.05


# =========================================================================
# Time utilities
# =========================================================================

class TestTime:
    def test_j2000_epoch(self):
        """J2000.0 = JD 2451545.0."""
        jd = julian_date(datetime(2000, 1, 1, 12, 0, 0))
        assert abs(jd - 2451545.0) < 1e-6

    def test_gmst_at_j2000(self):
        """GMST at J2000.0 epoch ~ 280.46 deg."""
        theta = np.degrees(gmst_rad(datetime(2000, 1, 1, 12, 0, 0)))
        assert abs(theta - 280.46) < 0.1


# =========================================================================
# Coordinate transforms
# =========================================================================

class TestCoordinates:
    def test_ecef_equator_prime_meridian(self):
        """(0, 0, 0) -> (a, 0, 0)."""
        r = geodetic_to_ecef(0.0, 0.0, 0.0)
        assert abs(r[0] - A_WGS84) < 1.0
        assert abs(r[1]) < 1.0
        assert abs(r[2]) < 1.0

    def test_ecef_north_pole(self):
        """(90, *, 0) -> (0, 0, b)."""
        r = geodetic_to_ecef(90.0, 0.0, 0.0)
        b = A_WGS84 * (1.0 - 1.0 / 298.257223563)
        assert abs(r[0]) < 1.0
        assert abs(r[1]) < 1.0
        assert abs(r[2] - b) < 1.0

    def test_enu_orthogonal(self):
        """ENU vectors are orthonormal."""
        e, n, u = enu_basis(48.14, 11.58)
        assert abs(np.dot(e, n)) < 1e-12
        assert abs(np.dot(e, u)) < 1e-12
        assert abs(np.dot(n, u)) < 1e-12
        assert abs(np.linalg.norm(e) - 1) < 1e-12
        assert abs(np.linalg.norm(n) - 1) < 1e-12
        assert abs(np.linalg.norm(u) - 1) < 1e-12

    def test_enu_right_handed(self):
        """E x N = U."""
        e, n, u = enu_basis(48.14, 11.58)
        cross = np.cross(e, n)
        np.testing.assert_allclose(cross, u, atol=1e-12)

    def test_vertical_axis(self):
        """zenith=0 gives the up vector."""
        e, n, u = enu_basis(48.14, 11.58)
        ax = measurement_axis(48.14, 11.58, 0.0, 0.0)
        np.testing.assert_allclose(ax, u, atol=1e-12)

    def test_horizontal_north(self):
        """zenith=90, azimuth=0 gives the north vector."""
        e, n, u = enu_basis(48.14, 11.58)
        ax = measurement_axis(48.14, 11.58, 90.0, 0.0)
        np.testing.assert_allclose(ax, n, atol=1e-12)

    def test_horizontal_east(self):
        """zenith=90, azimuth=90 gives the east vector."""
        e, n, u = enu_basis(48.14, 11.58)
        ax = measurement_axis(48.14, 11.58, 90.0, 90.0)
        np.testing.assert_allclose(ax, e, atol=1e-12)


# =========================================================================
# Ephemeris
# =========================================================================

class TestEphemeris:
    def test_sun_distance(self):
        """Sun distance ~ 1 AU +/- 3%."""
        dt = datetime(2025, 3, 20, 12, 0, 0)
        R = np.linalg.norm(sun_position_ecef(dt))
        assert abs(R / 1.496e11 - 1.0) < 0.03

    def test_moon_distance(self):
        """Moon distance ~ 385000 km +/- 7%."""
        dt = datetime(2025, 3, 20, 12, 0, 0)
        R = np.linalg.norm(moon_position_ecef(dt))
        assert abs(R / 3.844e8 - 1.0) < 0.07

    def test_moon_distance_varies(self):
        """Moon distance varies over a month (perigee/apogee)."""
        dists = []
        for day in range(30):
            dt = datetime(2025, 1, 1) + timedelta(days=day)
            dists.append(np.linalg.norm(moon_position_ecef(dt)))
        dists = np.array(dists)
        # Perigee-apogee range: ~40000-50000 km
        assert np.ptp(dists) > 3e7  # > 30000 km

    def test_sun_equinox_declination(self):
        """Near equinox, Sun declination ~ 0."""
        dt = datetime(2025, 3, 20, 12, 0, 0)
        pos = sun_position_ecef(dt)
        # In ECEF the z-component reflects declination;
        # at equinox |z/r| should be small
        r = np.linalg.norm(pos)
        dec = np.degrees(np.arcsin(pos[2] / r))
        assert abs(dec) < 2.0  # within 2 degrees of equator


# =========================================================================
# Tidal acceleration
# =========================================================================

class TestTidalAcceleration:
    def test_order_of_magnitude(self):
        """Vertical tidal acceleration ~ 1 um/s^2 for Moon."""
        dt = datetime(2025, 3, 20, 12, 0, 0)
        r = geodetic_to_ecef(48.14, 11.58, 500.0)
        R = moon_position_ecef(dt)
        a = tidal_acceleration(r, R, GM_MOON)
        mag = np.linalg.norm(a)
        assert 0.3e-6 < mag < 2.0e-6

    def test_sun_moon_ratio(self):
        """Solar tidal ~ 46% of lunar tidal (theoretical: ~0.46)."""
        dt = datetime(2025, 3, 20, 12, 0, 0)
        r = geodetic_to_ecef(48.14, 11.58, 500.0)
        a_m = np.linalg.norm(tidal_acceleration(r, moon_position_ecef(dt), GM_MOON))
        a_s = np.linalg.norm(tidal_acceleration(r, sun_position_ecef(dt), GM_SUN))
        ratio = a_s / a_m
        assert 0.3 < ratio < 0.7

    def test_zero_at_source_center(self):
        """Tidal acceleration vanishes when observer is at Earth center."""
        R = np.array([3.844e8, 0.0, 0.0])
        r = np.array([0.0, 0.0, 0.0])
        a = tidal_acceleration(r, R, GM_MOON)
        np.testing.assert_allclose(a, [0.0, 0.0, 0.0], atol=1e-20)


# =========================================================================
# Full pipeline: compute_g
# =========================================================================

class TestComputeG:
    def test_vertical_sensor(self):
        """g_total ~ 9.81 m/s^2 for vertical sensor."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           48.14, 11.58, 500.0)
        assert abs(result.g_total - 9.81) < 0.02

    def test_tidal_order_of_magnitude(self):
        """Tidal perturbation ~ um/s^2 level."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           48.14, 11.58, 500.0)
        assert abs(result.g_tidal) < 3e-6
        assert abs(result.g_tidal) > 1e-8

    def test_horizontal_smaller_static(self):
        """Horizontal sensor: g_static ~ 0."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           48.14, 11.58, 500.0,
                           zenith_deg=90.0, azimuth_deg=0.0)
        assert abs(result.g_static) < 1e-6

    def test_decomposition_sums(self):
        """g_tidal = g_tidal_moon + g_tidal_sun."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           48.14, 11.58, 500.0)
        assert abs(result.g_tidal - result.g_tidal_moon
                   - result.g_tidal_sun) < 1e-15

    def test_total_is_sum(self):
        """g_total = g_static + g_tidal."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           48.14, 11.58, 500.0)
        assert abs(result.g_total - result.g_static
                   - result.g_tidal) < 1e-15

    def test_gravimetric_amplification(self):
        """Tidal signal should be amplified by delta ~ 1.16 vs rigid Earth."""
        dt = datetime(2025, 3, 20, 12, 0, 0)
        r = geodetic_to_ecef(48.14, 11.58, 500.0)
        e_up = enu_basis(48.14, 11.58)[2]

        a_moon_rigid = np.dot(tidal_acceleration(r, moon_position_ecef(dt), GM_MOON), e_up)
        result = compute_g(dt, 48.14, 11.58, 500.0)

        # g_tidal_moon should be ~ delta * a_moon_rigid
        ratio = result.g_tidal_moon / a_moon_rigid
        assert abs(ratio - DELTA_GRAV) < 0.001


# =========================================================================
# Timeseries
# =========================================================================

class TestTimeseries:
    def test_length(self):
        """Correct number of time steps."""
        data = compute_timeseries(
            datetime(2025, 3, 20), datetime(2025, 3, 21),
            48.14, 11.58, 500.0, interval_minutes=60.0,
        )
        assert len(data.times) == 25  # 0, 1, ..., 24 hours
        assert len(data.g_total) == 25

    def test_tidal_variation(self):
        """Tidal signal varies over 48 hours (not constant)."""
        data = compute_timeseries(
            datetime(2025, 3, 20), datetime(2025, 3, 22),
            48.14, 11.58, 500.0, interval_minutes=30.0,
        )
        assert np.ptp(data.g_tidal) > 0.5e-6  # > 0.5 um/s^2

    def test_tidal_periods(self):
        """Strong spectral peaks near 12 h (semidiurnal) and/or 24 h (diurnal)."""
        data = compute_timeseries(
            datetime(2025, 3, 20), datetime(2025, 3, 27),
            48.14, 11.58, 500.0, interval_minutes=10.0,
        )
        signal = data.g_tidal - np.mean(data.g_tidal)
        n = len(signal)
        dt_sec = 600.0

        freqs = np.fft.rfftfreq(n, d=dt_sec)
        power = np.abs(np.fft.rfft(signal)) ** 2

        # Find power in semidiurnal band (11-13 h) and diurnal band (23-26 h)
        periods_h = np.where(freqs > 0, 1.0 / freqs / 3600.0, np.inf)
        sd_mask = (periods_h > 11.0) & (periods_h < 13.0)
        di_mask = (periods_h > 23.0) & (periods_h < 26.0)
        total_power = np.sum(power[1:])

        # Semidiurnal + diurnal should contain most of the tidal power
        tidal_power = np.sum(power[sd_mask]) + np.sum(power[di_mask])
        assert tidal_power / total_power > 0.5

    def test_static_constant(self):
        """Static gravity is constant across timeseries."""
        data = compute_timeseries(
            datetime(2025, 3, 20), datetime(2025, 3, 21),
            48.14, 11.58, 500.0,
        )
        assert np.ptp(data.g_static) == 0.0


# =========================================================================
# Edge cases
# =========================================================================

class TestEdgeCases:
    def test_south_pole(self):
        """Works at the south pole."""
        result = compute_g(datetime(2025, 6, 21, 12, 0, 0),
                           -90.0, 0.0, 2835.0)
        assert abs(result.g_total - 9.82) < 0.02

    def test_date_line(self):
        """Works at the date line (lon=180)."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           0.0, 180.0, 0.0)
        assert abs(result.g_total - 9.78) < 0.01

    def test_high_altitude(self):
        """Works at high altitude (airplane ~10 km)."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           48.14, 11.58, 10000.0)
        # Free-air: ~3 mGal/m -> 30 m/s^2 reduction, so g ~ 9.78
        assert 9.75 < result.g_total < 9.82

    def test_negative_longitude(self):
        """Works with negative longitude (western hemisphere)."""
        result = compute_g(datetime(2025, 3, 20, 12, 0, 0),
                           40.0, -74.0, 0.0)  # New York area
        assert abs(result.g_total - 9.80) < 0.02


# =========================================================================
# LabFrame and GravityField
# =========================================================================

class TestLabFrame:
    """Tests for LabFrame and GravityField."""

    DT = datetime(2025, 3, 20, 12, 0, 0)
    LAT, LON, ALT = 48.14, 11.58, 500.0

    @pytest.fixture
    def lab(self):
        return LabFrame(self.LAT, self.LON, self.ALT)

    @pytest.fixture
    def field(self, lab):
        return lab.field(self.DT)

    # -- Rotation vector --

    def test_omega_magnitude(self):
        """Omega magnitude = Earth rotation rate at any latitude."""
        for lat in [0.0, 30.0, 45.0, 60.0, 90.0]:
            lab = LabFrame(lat, 0.0)
            f = lab.field(self.DT)
            assert abs(np.linalg.norm(f.omega) - OMEGA) < 1e-15

    def test_omega_equator_points_north(self):
        """At equator, omega is purely North."""
        lab = LabFrame(0.0, 0.0)
        f = lab.field(self.DT)
        assert abs(f.omega[0]) < 1e-20  # no East
        assert abs(f.omega[1] - OMEGA) < 1e-15  # North
        assert abs(f.omega[2]) < 1e-15  # no Up

    def test_omega_pole_points_up(self):
        """At north pole, omega is purely Up."""
        lab = LabFrame(90.0, 0.0)
        f = lab.field(self.DT)
        assert abs(f.omega[0]) < 1e-20  # no East
        assert abs(f.omega[1]) < 1e-15  # no North
        assert abs(f.omega[2] - OMEGA) < 1e-15  # Up

    # -- Tidal gradient tensor --

    def test_tidal_gradient_symmetric(self, lab):
        """Tidal gradient tensor is symmetric."""
        from pytheas._core import (
            _tidal_gradient_tensor, _ecef_to_enu_tensor,
        )
        R_moon = moon_position_ecef(self.DT)
        T_ecef = _tidal_gradient_tensor(lab._r, R_moon, GM_MOON)
        T_enu = _ecef_to_enu_tensor(T_ecef, lab._enu)
        np.testing.assert_allclose(T_enu, T_enu.T, atol=1e-30)

    def test_tidal_gradient_traceless(self, lab):
        """Tidal gradient tensor is traceless (Laplace equation)."""
        from pytheas._core import (
            _tidal_gradient_tensor, _ecef_to_enu_tensor,
        )
        R_moon = moon_position_ecef(self.DT)
        T_ecef = _tidal_gradient_tensor(lab._r, R_moon, GM_MOON)
        T_enu = _ecef_to_enu_tensor(T_ecef, lab._enu)
        assert abs(np.trace(T_enu)) < 1e-25

    def test_tidal_gradient_magnitude(self, lab):
        """Lunar tidal gradient ~ GM_Moon/d^3 ~ 1e-13 s^-2 (~0.1 milliEotvos)."""
        from pytheas._core import (
            _tidal_gradient_tensor, _ecef_to_enu_tensor,
        )
        R_moon = moon_position_ecef(self.DT)
        T_ecef = DELTA_GRAV * _tidal_gradient_tensor(
            lab._r, R_moon, GM_MOON)
        T_enu = _ecef_to_enu_tensor(T_ecef, lab._enu)
        # Max component ~ 2*GM/d^3 * delta ~ 2e-13 s^-2
        max_comp = np.max(np.abs(T_enu))
        assert 1e-14 < max_comp < 1e-12

    # -- Earth gradient tensor --

    def test_earth_gradient_trace(self, lab):
        """Earth gradient tensor trace = -2*Omega^2 (Poisson condition)."""
        from pytheas._core import _earth_gradient_tensor
        T_earth = _earth_gradient_tensor(self.LAT, self.ALT)
        expected = -2.0 * OMEGA ** 2
        assert abs(np.trace(T_earth) - expected) < 1e-20

    def test_earth_gradient_T_UU(self):
        """Free-air gradient T_UU ~ -3086 Eotvos at mid-latitude, sea level."""
        from pytheas._core import _earth_gradient_tensor
        T = _earth_gradient_tensor(45.0, 0.0)
        T_UU_eotvos = T[2, 2] * 1e9  # convert s^-2 to nanoEotvos... no
        T_UU_eotvos = T[2, 2] / 1e-9  # convert s^-2 to Eotvos
        # T_UU should be about -3086e-9 s^-2 = -3086 nE = -3.086 uGal/m
        # In Eotvos (1E = 1e-9 s^-2): ~ -3086 E
        assert abs(T_UU_eotvos - (-3086)) < 20

    def test_earth_gradient_symmetric(self):
        """Earth gradient tensor is symmetric."""
        from pytheas._core import _earth_gradient_tensor
        T = _earth_gradient_tensor(self.LAT, self.ALT)
        np.testing.assert_allclose(T, T.T, atol=1e-30)

    # -- Full tensor --

    def test_total_tensor_trace(self, field):
        """Total tensor trace = -2*Omega^2 (tidal is traceless)."""
        expected = -2.0 * OMEGA ** 2
        assert abs(np.trace(field.T) - expected) < 1e-20

    def test_total_tensor_symmetric(self, field):
        """Total gravity gradient tensor is symmetric."""
        np.testing.assert_allclose(field.T, field.T.T, atol=1e-25)

    # -- GravityField methods --

    def test_at_zero_returns_g(self, field):
        """at([0,0,0]) returns g."""
        np.testing.assert_allclose(field.at([0, 0, 0]), field.g, atol=1e-20)

    def test_at_nonzero(self, field):
        """at(dx) = g + T @ dx."""
        dx = np.array([1.0, 0.0, 0.0])
        expected = field.g + field.T @ dx
        np.testing.assert_allclose(field.at(dx), expected, atol=1e-20)

    def test_eom_zero_returns_g(self, field):
        """eom(0, 0) returns g."""
        np.testing.assert_allclose(
            field.eom([0, 0, 0], [0, 0, 0]), field.g, atol=1e-20)

    def test_coriolis_direction(self, field):
        """Coriolis deflects eastward-moving object downward at mid-lat."""
        # v = v_east, Coriolis = -2*(omega x v)
        # At mid-lat omega has N and U components.
        # omega x v_east = (omega_N * 0 - omega_U * 0, omega_U * v - 0,
        #                    0 - omega_N * v) -- using ENU cross product
        # Actually: omega = (0, omega_N, omega_U), v = (v, 0, 0)
        # omega x v = (omega_N*0 - omega_U*0, omega_U*v - 0*0,
        #              0*0 - omega_N*v) = (0, omega_U*v, -omega_N*v)
        # Coriolis = -2*(0, omega_U*v, -omega_N*v)
        #          = (0, -2*omega_U*v, 2*omega_N*v)
        # So Up component is positive (deflected upward? No, let me check)
        # Actually for eastward motion the vertical Coriolis is:
        # a_U = 2*omega_N*v > 0 (upward for eastward motion)
        v_east = 1.0
        a = field.eom([0, 0, 0], [v_east, 0, 0])
        a_no_v = field.eom([0, 0, 0], [0, 0, 0])
        coriolis = a - a_no_v
        # North component: -2*omega_U*v < 0 (southward deflection)
        assert coriolis[1] < 0
        # Up component: 2*omega_N*v > 0 (Eotvos effect)
        assert coriolis[2] > 0

    def test_reading_vertical(self, field):
        """Vertical reading ~ -g_normal."""
        r = field.reading([0, 0, 1])
        assert abs(r - (-field.g_normal)) < 1e-4

    # -- Shapes and types --

    def test_field_shapes(self, field):
        """All arrays have correct shapes."""
        assert field.g.shape == (3,)
        assert field.omega.shape == (3,)
        assert field.T.shape == (3, 3)
        assert field.g_tidal_moon.shape == (3,)
        assert field.g_tidal_sun.shape == (3,)
        assert isinstance(field.g_normal, float)

    # -- Consistency with compute_g --

    def test_consistency_with_compute_g(self, lab, field):
        """LabFrame tidal components match compute_g component-by-component."""
        result = compute_g(self.DT, self.LAT, self.LON, self.ALT)
        # Normal gravity must agree exactly
        assert abs(field.g_normal - result.g_normal) < 1e-15
        # Tidal components in Up direction must match (ENU[2] = ECEF dot e_up)
        assert abs(field.g_tidal_moon[2] - result.g_tidal_moon) < 1e-15
        assert abs(field.g_tidal_sun[2] - result.g_tidal_sun) < 1e-15
        # Cross-check: g_total = g_normal + tidal_up, g[2] = -g_normal + tidal_up
        # So g_total + g[2] = 2 * tidal_up
        tidal_up = field.g_tidal_moon[2] + field.g_tidal_sun[2]
        assert abs((result.g_total + field.g[2]) - 2 * tidal_up) < 1e-14

    # -- Timeseries --

    def test_timeseries_length(self, lab):
        """Timeseries returns correct number of fields."""
        ts = lab.timeseries(
            datetime(2025, 3, 20), datetime(2025, 3, 21),
            interval_minutes=60.0)
        assert len(ts['times']) == 25
        assert len(ts['fields']) == 25

    def test_timeseries_returns_gravity_fields(self, lab):
        """Timeseries returns GravityField objects."""
        ts = lab.timeseries(
            datetime(2025, 3, 20), datetime(2025, 3, 20, 1, 0),
            n_samples=3)
        for f in ts['fields']:
            assert isinstance(f, GravityField)

    def test_timeseries_n_samples(self, lab):
        """n_samples overrides interval_minutes."""
        ts = lab.timeseries(
            datetime(2025, 3, 20), datetime(2025, 3, 21),
            n_samples=10)
        assert len(ts['times']) == 10
        assert len(ts['fields']) == 10


# =========================================================================
# Dataclass behavior
# =========================================================================

class TestDataclassAPI:
    """Tests for GravityResult and TimeSeries dataclass behavior."""

    DT = datetime(2025, 3, 20, 12, 0, 0)

    def test_compute_g_returns_gravity_result(self):
        result = compute_g(self.DT, 48.14, 11.58, 500.0)
        assert isinstance(result, GravityResult)

    def test_compute_timeseries_returns_timeseries(self):
        data = compute_timeseries(
            self.DT, self.DT + timedelta(hours=1),
            48.14, 11.58, 500.0, n_samples=3)
        assert isinstance(data, TimeSeries)

    def test_gravity_result_frozen(self):
        result = compute_g(self.DT, 48.14, 11.58, 500.0)
        with pytest.raises(FrozenInstanceError):
            result.g_total = 0.0

    def test_timeseries_frozen(self):
        data = compute_timeseries(
            self.DT, self.DT + timedelta(hours=1),
            48.14, 11.58, 500.0, n_samples=3)
        with pytest.raises(FrozenInstanceError):
            data.g_normal = 0.0

    def test_gravity_result_asdict(self):
        result = compute_g(self.DT, 48.14, 11.58, 500.0)
        d = asdict(result)
        assert isinstance(d, dict)
        assert set(d.keys()) == {
            'g_total', 'g_static', 'g_normal', 'cos_zenith',
            'g_tidal', 'g_tidal_moon', 'g_tidal_sun',
        }
        assert d['g_total'] == result.g_total

    def test_timeseries_asdict(self):
        data = compute_timeseries(
            self.DT, self.DT + timedelta(hours=1),
            48.14, 11.58, 500.0, n_samples=3)
        d = asdict(data)
        assert isinstance(d, dict)
        assert 'times' in d
        assert 'g_total' in d
        assert d['g_normal'] == data.g_normal


# =========================================================================
# External validation (1b)
# =========================================================================

class TestExternalValidation:
    """Validate against JPL Horizons (DE441) and GRS80 reference values.

    Reference positions: JPL Horizons API, geocentric ICRF (J2000 equatorial).
    Since ECI→ECEF only rotates around z, the z-component and geocentric
    distance are frame-independent checks.
    """

    # JPL Horizons DE441 Moon geocentric ICRF positions (km)
    MOON_REFS = [
        # (datetime, distance_km, z_km)
        (datetime(2024, 1, 1, 0, 0, 0),   404667.519, 89342.281),
        (datetime(2024, 6, 21, 12, 0, 0),  382270.712, -179591.076),
        (datetime(2025, 3, 20, 12, 0, 0),  401799.262, -179121.060),
        (datetime(2025, 9, 23, 0, 0, 0),   400463.745, -50275.994),
    ]

    # JPL Horizons DE441 Sun geocentric ICRF positions (km)
    SUN_REFS = [
        # (datetime, distance_km, z_km)
        (datetime(2024, 1, 1, 0, 0, 0),   147102329.179, -57668106.240),
        (datetime(2024, 6, 21, 12, 0, 0),  152026590.092, 60463868.434),
        (datetime(2025, 3, 20, 12, 0, 0),  148988141.499, -233235.764),
        (datetime(2025, 9, 23, 0, 0, 0),   150130470.991, 129666.881),
    ]

    # GRS80/WGS84 normal gravity on the ellipsoid (h=0)
    GRAVITY_REFS = [
        # (lat_deg, gamma_m_s2)
        (0.0,  9.7803253141),   # equatorial (exact GRS80)
        (45.0, 9.8061992032),   # mid-latitude (Somigliana)
        (90.0, 9.8321849378),   # polar (exact GRS80)
    ]

    @pytest.mark.parametrize("dt,ref_dist_km,ref_z_km", MOON_REFS,
                             ids=["2024-01-01", "2024-06-21", "2025-03-20",
                                  "2025-09-23"])
    def test_moon_distance(self, dt, ref_dist_km, ref_z_km):
        """Moon geocentric distance within 200 km of DE441."""
        pos = moon_position_ecef(dt)
        dist_km = np.linalg.norm(pos) / 1e3
        assert abs(dist_km - ref_dist_km) < 200, (
            f"Moon distance: {dist_km:.1f} km vs {ref_dist_km:.1f} km "
            f"(err = {abs(dist_km - ref_dist_km):.1f} km)")

    @pytest.mark.parametrize("dt,ref_dist_km,ref_z_km", MOON_REFS,
                             ids=["2024-01-01", "2024-06-21", "2025-03-20",
                                  "2025-09-23"])
    def test_moon_direction(self, dt, ref_dist_km, ref_z_km):
        """Moon declination within 0.15 deg of DE441."""
        pos = moon_position_ecef(dt)
        # z_ecef = z_eci (rotation about z-axis preserves z)
        sin_dec = pos[2] / np.linalg.norm(pos) / 1e3 * 1e3
        sin_dec_ref = ref_z_km / ref_dist_km
        # Angular difference in degrees
        dec = np.degrees(np.arcsin(np.clip(sin_dec, -1, 1)))
        dec_ref = np.degrees(np.arcsin(np.clip(sin_dec_ref, -1, 1)))
        assert abs(dec - dec_ref) < 0.15, (
            f"Moon declination: {dec:.3f} deg vs {dec_ref:.3f} deg")

    @pytest.mark.parametrize("dt,ref_dist_km,ref_z_km", SUN_REFS,
                             ids=["2024-01-01", "2024-06-21", "2025-03-20",
                                  "2025-09-23"])
    def test_sun_distance(self, dt, ref_dist_km, ref_z_km):
        """Sun geocentric distance within 10000 km (~0.007%) of DE441."""
        pos = sun_position_ecef(dt)
        dist_km = np.linalg.norm(pos) / 1e3
        assert abs(dist_km - ref_dist_km) < 10000, (
            f"Sun distance: {dist_km:.1f} km vs {ref_dist_km:.1f} km "
            f"(err = {abs(dist_km - ref_dist_km):.1f} km)")

    @pytest.mark.parametrize("dt,ref_dist_km,ref_z_km", SUN_REFS,
                             ids=["2024-01-01", "2024-06-21", "2025-03-20",
                                  "2025-09-23"])
    def test_sun_direction(self, dt, ref_dist_km, ref_z_km):
        """Sun declination within 0.2 deg of DE441.

        Tolerance accounts for Meeus low-precision (~1 arcmin) ephemeris,
        especially near equinoxes where z/|r| -> 0.
        """
        pos = sun_position_ecef(dt)
        sin_dec = pos[2] / np.linalg.norm(pos)
        sin_dec_ref = ref_z_km / ref_dist_km
        dec = np.degrees(np.arcsin(np.clip(sin_dec, -1, 1)))
        dec_ref = np.degrees(np.arcsin(np.clip(sin_dec_ref, -1, 1)))
        assert abs(dec - dec_ref) < 0.2, (
            f"Sun declination: {dec:.4f} deg vs {dec_ref:.4f} deg")

    @pytest.mark.parametrize("lat_deg,ref_gamma", GRAVITY_REFS,
                             ids=["equator", "45deg", "pole"])
    def test_normal_gravity(self, lat_deg, ref_gamma):
        """Normal gravity matches GRS80 reference to < 1 mGal."""
        gamma = normal_gravity(lat_deg, 0.0)
        assert abs(gamma - ref_gamma) < 1e-5, (
            f"gamma({lat_deg}) = {gamma:.10f} vs {ref_gamma:.10f} "
            f"(err = {abs(gamma - ref_gamma) * 1e6:.2f} uGal)")


# =========================================================================
# Physics properties (1c)
# =========================================================================

class TestPhysicsProperties:
    """Property-based tests verifying physical invariants."""

    DT = datetime(2025, 3, 20, 12, 0, 0)

    def test_tidal_vanishes_at_large_distance(self):
        """Tidal acceleration -> 0 as body distance -> infinity."""
        r = geodetic_to_ecef(48.14, 11.58, 500.0)
        R_moon = moon_position_ecef(self.DT)
        R_dir = R_moon / np.linalg.norm(R_moon)

        a_1x = np.linalg.norm(tidal_acceleration(r, R_moon, GM_MOON))
        a_10x = np.linalg.norm(
            tidal_acceleration(r, 10 * R_moon, GM_MOON))
        a_100x = np.linalg.norm(
            tidal_acceleration(r, 100 * R_moon, GM_MOON))

        # Tidal ~ 1/R^3, so at 10x distance should be ~1000x smaller
        assert a_10x < a_1x / 500
        assert a_100x < a_1x / 500000

    def test_tidal_gradient_limit(self):
        """Near Earth center, tidal approaches gradient approximation.

        For r/R << 1: a_tidal ~ -GM/R^3 * r + 3*GM*(R_hat . r)*R_hat/R^3
        """
        R = np.array([3.844e8, 0.0, 0.0])  # Moon-like distance on x-axis
        R_norm = np.linalg.norm(R)
        R_hat = R / R_norm

        # Small displacement from origin
        r_small = np.array([1000.0, 500.0, 200.0])  # 1 km offset

        a_exact = tidal_acceleration(r_small, R, GM_MOON)

        # Gradient approximation
        a_grad = GM_MOON / R_norm**3 * (
            -r_small + 3.0 * np.dot(R_hat, r_small) * R_hat)

        # Should agree to high precision (r/R ~ 3e-6)
        np.testing.assert_allclose(a_exact, a_grad, rtol=1e-5)

    @pytest.mark.parametrize("lat,lon", [
        (0.0, 0.0),       # equator, prime meridian
        (90.0, 0.0),      # north pole
        (-90.0, 0.0),     # south pole
        (0.0, 180.0),     # equator, date line
        (0.0, -90.0),     # equator, 90W
        (-33.87, 151.21), # Sydney
        (35.68, 139.69),  # Tokyo
    ])
    def test_enu_right_handed(self, lat, lon):
        """ENU basis is right-handed (E x N = U) at all locations."""
        e, n, u = enu_basis(lat, lon)
        cross = np.cross(e, n)
        np.testing.assert_allclose(cross, u, atol=1e-12)
        # Also orthonormal
        assert abs(np.linalg.norm(e) - 1) < 1e-12
        assert abs(np.linalg.norm(n) - 1) < 1e-12
        assert abs(np.linalg.norm(u) - 1) < 1e-12

    def test_normal_gravity_bounds(self):
        """gamma_e <= gamma(lat, 0) <= gamma_p for all latitudes."""
        rng = np.random.default_rng(42)
        lats = rng.uniform(-90, 90, 50)
        for lat in lats:
            g = normal_gravity(lat, 0.0)
            assert GAMMA_E <= g <= GAMMA_P, (
                f"gamma({lat:.2f}) = {g:.10f} outside "
                f"[{GAMMA_E}, {GAMMA_P}]")

    @pytest.mark.parametrize("lat", [0.0, 30.0, 45.0, 60.0, 90.0, -45.0])
    def test_gravity_gradient_symmetric(self, lat):
        """Total gravity gradient tensor is symmetric at any latitude."""
        from pytheas._core import _earth_gradient_tensor
        T = _earth_gradient_tensor(lat, 0.0)
        np.testing.assert_allclose(T, T.T, atol=1e-30)

    def test_coriolis_perpendicular_to_velocity(self):
        """Coriolis acceleration is always perpendicular to velocity."""
        lab = LabFrame(48.14, 11.58, 500.0)
        field = lab.field(self.DT)

        rng = np.random.default_rng(123)
        for _ in range(20):
            v = rng.standard_normal(3) * 10.0  # random velocity, ~10 m/s
            a_total = field.eom([0, 0, 0], v)
            a_no_v = field.eom([0, 0, 0], [0, 0, 0])
            coriolis = a_total - a_no_v  # = -2*(omega x v)
            # Coriolis must be perpendicular to v (within floating-point)
            dot = np.dot(coriolis, v)
            assert abs(dot) < 1e-12, (
                f"Coriolis not perpendicular to v: dot = {dot}")
