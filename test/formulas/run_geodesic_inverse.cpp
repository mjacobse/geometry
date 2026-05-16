#include <geometry_test_common.hpp>

#include <boost/geometry/formulas/andoyer_inverse.hpp>
#include <boost/geometry/srs/srs.hpp>

#include <fstream>
#include <iomanip>

int test_main(int, char*[])
{
    const boost::geometry::srs::spheroid<double> boost_wgs84{};
    using BoostInverse = boost::geometry::formula::andoyer_inverse<double, true, true, true, true, true>;

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

    std::ifstream test_data("GeodTest.dat");
    double lat1;
    double lon1;
    double azi1;
    double lat2;
    double lon2;
    double azi2;
    double s12;
    double sigma12;
    double m12;
    double S12;
    while( test_data >> lat1 >> lon1 >> azi1 >> lat2 >> lon2 >> azi2 >>
           s12 >> sigma12 >> m12 >> S12 ) {
        double const d2r = boost::geometry::math::d2r<double>();
        double const r2d = boost::geometry::math::r2d<double>();
        boost::geometry::formula::result_inverse result = BoostInverse::apply(lon1 * d2r, lat1 * d2r, lon2 * d2r, lat2 * d2r, boost_wgs84);

        results
            << lat1 << ',' << lon1 << ','
            << lat2 << ',' << lon2 << ','
            << result.distance << ',' << s12 << ','
            << result.azimuth * r2d << ',' << azi1 << ','
            << result.reverse_azimuth * r2d << ',' << azi2 << ','
            << result.reduced_length << ',' << m12 << ','
            << result.geodesic_scale << ',' << result.geodesic_scale << ','  // not in test data
            << '\n';
    }

    return 0;
}
