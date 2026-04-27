#include <geometry_test_common.hpp>

#include <boost/geometry/formulas/andoyer_inverse.hpp>
#include <boost/geometry/srs/srs.hpp>

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>

#include <fstream>
#include <iomanip>
#include <random>

int test_main(int, char*[])
{
    const GeographicLib::Geodesic geographiclib_wgs84 = GeographicLib::Geodesic::WGS84();

    const boost::geometry::srs::spheroid<double> boost_wgs84{};
    using BoostInverse = boost::geometry::formula::andoyer_inverse<double, true, true, true, true, true>;

    std::mt19937 rng{0};

    std::ofstream results("results.csv");
    results << std::setprecision(17);
    results
        << "lat1,lon1,lat2,lon2,"
        << "boost_distance,geographiclib_distance,"
        << "boost_azimuth,geographiclib_azimuth,"
        << "boost_reverse_azimuth,geographiclib_reverse_azimuth,"
        << "boost_reduced_length,geographiclib_reduced_length,"
        << "boost_geodesic_scale,geographiclib_geodesic_scale,"
        << "\n";
    for (int i = 0; i < 100000; ++i)
    {
        double const lat1 = std::uniform_real_distribution{-90.0, 90.0}(rng);
        double const lon1 = std::uniform_real_distribution{-180.0, 180.0}(rng);
        double const sigma_lat = 90.0 * std::pow(10.0, std::uniform_real_distribution{-16.0, 0.0}(rng));
        double const sigma_lon = 180.0 * std::pow(10.0, std::uniform_real_distribution{-16.0, 0.0}(rng));
        double lat2;
        double lon2;
        do { lat2 = lat1 + std::normal_distribution{0.0, sigma_lat}(rng); } while (std::abs(lat2) > 90.0);
        do { lon2 = lon1 + std::normal_distribution{0.0, sigma_lon}(rng); } while (std::abs(lon2) > 180.0);

        double const d2r = boost::geometry::math::d2r<double>();
        double const r2d = boost::geometry::math::r2d<double>();

        boost::geometry::formula::result_inverse result = BoostInverse::apply(lon1 * d2r, lat1 * d2r, lon2 * d2r, lat2 * d2r, boost_wgs84);

        double s12, azi1, azi2, m12, M12, M21, S12;
        geographiclib_wgs84.Inverse(lat1, lon1, lat2, lon2, s12, azi1, azi2, m12, M12, M21, S12);

        results
            << lat1 << ',' << lon1 << ','
            << lat2 << ',' << lon2 << ','
            << result.distance << ',' << s12 << ','
            << result.azimuth * r2d << ',' << azi1 << ','
            << result.reverse_azimuth * r2d << ',' << azi2 << ','
            << result.reduced_length << ',' << m12 << ','
            << result.geodesic_scale << ',' << M12 << ','
            << '\n';
    }

    return 0;
}
