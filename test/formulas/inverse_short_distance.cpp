// Short-distance accuracy test for andoyer_inverse.
// Reference for the regression: https://github.com/boostorg/geometry/issues/1217
// Fix proposed in PR #1461.
//
// On develop the issue_1217 case fails: andoyer returns ~9 cm for a true
// distance of ~0.67 mm. After the PR, andoyer agrees with vincenty to
// well under a millimetre across the whole range tested here.

#include <geometry_test_common.hpp>

#include <boost/geometry/formulas/andoyer_inverse.hpp>
#include <boost/geometry/formulas/vincenty_inverse.hpp>
#include <boost/geometry/srs/spheroid.hpp>
#include <boost/geometry/util/math.hpp>

namespace
{

// Maximum allowed disagreement between andoyer and vincenty for the cases
// below. 1 mm comfortably catches the ~9 cm regression from issue #1217 and
// the cancellation-to-zero failures, while leaving headroom for the
// first-order Andoyer-Lambert approximation error at the metre scale.
double const tolerance_m = 1.0e-3;

struct short_case
{
    char const* name;
    double lon1, lat1, lon2, lat2;
    double distance_expected;
};

void test_short(short_case const& c, double abs_tol_meters = tolerance_m)
{
    double const d2r = bg::math::d2r<double>();
    bg::srs::spheroid<double> const spheroid; // WGS84 by default

    using andoyer_t = bg::formula::andoyer_inverse<double, true, false, false, false, false>;
    //using thomas_t = bg::formula::thomas_inverse<double, true, false, false, false, false>;
    using vincenty_t = bg::formula::vincenty_inverse<double, true, false, false, false, false>;
    using karney_t = bg::formula::karney_inverse<double, true, false, false, false, false>;

    double const lon1r = c.lon1 * d2r;
    double const lat1r = c.lat1 * d2r;
    double const lon2r = c.lon2 * d2r;
    double const lat2r = c.lat2 * d2r;

    double const distance_andoyer = andoyer_t::apply(lon1r, lat1r, lon2r, lat2r, spheroid).distance;
    //double const distance_thomas = thomas_t::apply(lon1r, lat1r, lon2r, lat2r, spheroid).distance;
    double const distance_vincenty = vincenty_t::apply(lon1r, lat1r, lon2r, lat2r, spheroid).distance;
    double const distance_karney = karney_t::apply(lon1r, lat1r, lon2r, lat2r, spheroid).distance;

    double const percent_tolerance = 1.0;  // allow error of 1%
    BOOST_TEST_INFO_SCOPE(c.name);
    BOOST_CHECK_CLOSE(distance_andoyer, c.distance_expected, 1e-8);
    //BOOST_CHECK_CLOSE(distance_thomas, c.distance_expected, percent_tolerance); // TODO: Thomas is very inaccurate
    BOOST_CHECK_CLOSE(distance_vincenty, c.distance_expected, 1e-8);
    BOOST_CHECK_CLOSE(distance_karney, c.distance_expected, percent_tolerance);
}

} // namespace

int test_main(int, char*[])
{
    // reference expected distance values obtained with GeodSolve/GeographicLib

    // Marquee case from issue #1217: ~0.67 mm true distance at mid-latitude.
    // Develop returns ~9 cm here (the regression the PR fixes).
    test_short({"issue_1217",      8.81, 53.08, 8.81000001, 53.08, 0.0006701306 });

    // East-west steps at the equator, sweeping sub-mm to ~10 m.
    test_short({"sub_mm_equator",  0.0,   0.0,  1.0e-8,      0.0, 0.0011131949 });
    test_short({"cm_equator",      0.0,   0.0,  1.0e-7,      0.0, 0.0111319491 });
    test_short({"m_equator",       0.0,   0.0,  1.0e-5,      0.0, 1.1131949079 });
    test_short({"10m_equator",     0.0,   0.0,  1.0e-4,      0.0, 11.1319490793 });

    // North-south steps along a meridian.
    test_short({"sub_mm_meridian", 0.0,  45.0,  0.0,         45.0 + 1.0e-8, 0.0011113192 });
    test_short({"m_meridian",      0.0,  45.0,  0.0,         45.0 + 1.0e-5, 1.1113177758 });

    // High-latitude oblique step (cos(lat) is small, so longitude differences shrink).
    test_short({"oblique_70N",     10.0, 70.0, 10.0 + 1.0e-7, 70.0 + 1.0e-7, 0.0117916474 });

    return 0;
}
