// Boost.Geometry - gis-projections (based on PROJ4)

// Copyright (c) 2008-2015 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2017, 2018, 2019.
// Modifications copyright (c) 2017-2019, Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// This file is converted from PROJ4, http://trac.osgeo.org/proj
// PROJ4 is originally written by Gerald Evenden (then of the USGS)
// PROJ4 is maintained by Frank Warmerdam
// PROJ4 is converted to Boost.Geometry by Barend Gehrels

// Last updated version of proj: 5.0.0

// Original copyright notice:

// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#ifndef BOOST_GEOMETRY_PROJECTIONS_ORTHO_HPP
#define BOOST_GEOMETRY_PROJECTIONS_ORTHO_HPP

#include <boost/config.hpp>
#include <boost/geometry/util/math.hpp>
#include <boost/math/special_functions/hypot.hpp>

#include <boost/geometry/srs/projections/impl/base_static.hpp>
#include <boost/geometry/srs/projections/impl/base_dynamic.hpp>
#include <boost/geometry/srs/projections/impl/projects.hpp>
#include <boost/geometry/srs/projections/impl/factory_entry.hpp>

namespace boost { namespace geometry
{

namespace projections
{
    #ifndef DOXYGEN_NO_DETAIL
    namespace detail { namespace ortho
    {

            enum mode_type {
                n_pole = 0,
                s_pole = 1,
                equit  = 2,
                obliq  = 3
            };

            template <typename T>
            struct par_ortho
            {
                T   sinph0;
                T   cosph0;
                T   nu0;
                mode_type mode;
            };

            static const double epsilon10 = 1.e-10;

            template <typename T, typename Parameters>
            struct base_ortho_spheroid
            {
                par_ortho<T> m_proj_parm;

                // FORWARD(s_forward)  spheroid
                // Project coordinates from geographic (lon, lat) to cartesian (x, y)
                inline void fwd(Parameters const& par, T const& lp_lon, T const& lp_lat, T& xy_x, T& xy_y) const
                {
                    static const T half_pi = detail::half_pi<T>();

                    T coslam, cosphi, sinphi;

                    cosphi = cos(lp_lat);
                    coslam = cos(lp_lon);
                    switch (this->m_proj_parm.mode) {
                    case equit:
                        if (cosphi * coslam < - epsilon10) {
                            BOOST_THROW_EXCEPTION( projection_exception(error_tolerance_condition) );
                        }
                        xy_y = sin(lp_lat);
                        break;
                    case obliq:
                        if (this->m_proj_parm.sinph0 * (sinphi = sin(lp_lat)) +
                           this->m_proj_parm.cosph0 * cosphi * coslam < - epsilon10) {
                            BOOST_THROW_EXCEPTION( projection_exception(error_tolerance_condition) );
                        }
                        xy_y = this->m_proj_parm.cosph0 * sinphi - this->m_proj_parm.sinph0 * cosphi * coslam;
                        break;
                    case n_pole:
                        coslam = - coslam;
                        BOOST_FALLTHROUGH;
                    case s_pole:
                        if (fabs(lp_lat - par.phi0) - epsilon10 > half_pi) {
                            BOOST_THROW_EXCEPTION( projection_exception(error_tolerance_condition) );
                        }
                        xy_y = cosphi * coslam;
                        break;
                    }
                    xy_x = cosphi * sin(lp_lon);
                }

                // INVERSE(s_inverse)  spheroid
                // Project coordinates from cartesian (x, y) to geographic (lon, lat)
                inline void inv(Parameters const& par, T xy_x, T xy_y, T& lp_lon, T& lp_lat) const
                {
                    static const T half_pi = detail::half_pi<T>();

                    T rh, cosc, sinc;

                    if ((sinc = (rh = boost::math::hypot(xy_x, xy_y))) > 1.) {
                        if ((sinc - 1.) > epsilon10) {
                            BOOST_THROW_EXCEPTION( projection_exception(error_tolerance_condition) );
                        }
                        sinc = 1.;
                    }
                    cosc = sqrt(1. - sinc * sinc); /* in this range OK */
                    if (fabs(rh) <= epsilon10) {
                        lp_lat = par.phi0;
                        lp_lon = 0.0;
                    } else {
                        switch (this->m_proj_parm.mode) {
                        case n_pole:
                            xy_y = -xy_y;
                            lp_lat = acos(sinc);
                            break;
                        case s_pole:
                            lp_lat = - acos(sinc);
                            break;
                        case equit:
                            lp_lat = xy_y * sinc / rh;
                            xy_x *= sinc;
                            xy_y = cosc * rh;
                            goto sinchk;
                        case obliq:
                            lp_lat = cosc * this->m_proj_parm.sinph0 + xy_y * sinc * this->m_proj_parm.cosph0 /rh;
                            xy_y = (cosc - this->m_proj_parm.sinph0 * lp_lat) * rh;
                            xy_x *= sinc * this->m_proj_parm.cosph0;
                        sinchk:
                            if (fabs(lp_lat) >= 1.)
                                lp_lat = lp_lat < 0. ? -half_pi : half_pi;
                            else
                                lp_lat = asin(lp_lat);
                            break;
                        }
                        lp_lon = (xy_y == 0. && (this->m_proj_parm.mode == obliq || this->m_proj_parm.mode == equit))
                             ? (xy_x == 0. ? 0. : xy_x < 0. ? -half_pi : half_pi)
                                           : atan2(xy_x, xy_y);
                    }
                }

                static inline std::string get_name()
                {
                    return "ortho_spheroid";
                }

            };

            template <typename T, typename Parameters>
            struct base_ortho_ellipsoid
            {
                par_ortho<T> m_proj_parm;

                // FORWARD(s_forward)  ellipsoid
                // Project coordinates from geographic (lon, lat) to cartesian (x, y)
                inline void fwd(Parameters const &par, T const &lp_lon, T const &lp_lat, T &xy_x, T &xy_y) const
                {
                    const T cosphi = cos(lp_lat);
                    const T sinphi = sin(lp_lat);
                    const T coslam = cos(lp_lon);
                    const T sinlam = sin(lp_lon);
                    const T nu = 1.0 / sqrt(1.0 - par.es * sinphi * sinphi);
                    switch (this->m_proj_parm.mode)
                    {
                    case equit:
                        // TODO: out of bounds check
                        xy_y = (1.0 - par.es) * nu * sinphi;
                        break;
                    case obliq:
                        {
                        // TODO: out of bounds check
                        const T cosphi0 = this->m_proj_parm.cosph0;
                        const T sinphi0 = this->m_proj_parm.sinph0;
                        const T nu0 = this->m_proj_parm.nu0;
                        xy_y = nu * (sinphi * cosphi0 - cosphi * sinphi0 * coslam) + par.es * (nu0 * sinphi0 - nu * sinphi) * cosphi0;
                        break;
                        }
                    case n_pole:
                        // TODO: out of bounds check
                        xy_y = -nu * cosphi * coslam;
                        break;
                    case s_pole:
                        // TODO: out of bounds check
                        xy_y = nu * cosphi * coslam;
                        break;
                    }
                    xy_x = nu * cosphi * sinlam;
                }

                // INVERSE(s_inverse)  ellipsoid
                // Project coordinates from cartesian (x, y) to geographic (lon, lat)
                inline void inv(Parameters const& par, T xy_x, T xy_y, T& lp_lon, T& lp_lat) const
                {
                    static const T half_pi = detail::half_pi<T>();

                    if (this->m_proj_parm.mode == n_pole || this->m_proj_parm.mode == s_pole) {
                        // x = sign * nu * cosphi * sinlam
                        // y = nu * cosphi * coslam
                        // (1)/(2) => sign * x/y = tan(lam) => lam = atan2(sign * x, y)
                        // (1)^2 + (2)^2
                        // => x^2 + y^2 = nu^2 * cosphi^2
                        // => (x^2 + y^2) * (1 - e^2 * sinphi^2) = cosphi^2
                        // => (x^2 + y^2) * (1 - e^2 * (1 - cosphi^2)) = cosphi^2
                        // => (x^2 + y^2) * (1 - e^2) = cosphi^2 (1 - e^2 (x^2 + y^2))
                        // => (x^2 + y^2) * (1 - e^2) / (1 - e^2 (x^2 + y^2)) = cosphi^2

                        const T rh_squared = boost::math::pow<2>(xy_x) + boost::math::pow<2>(xy_y);
                        if (rh_squared > 1.0) {
                            BOOST_THROW_EXCEPTION( projection_exception(error_tolerance_condition) );
                        }

                        const T sign = this->m_proj_parm.mode == n_pole ? -1 : 1;
                        lp_lon = atan2(xy_x, sign * xy_y);
                        lp_lat = acos(sqrt(rh_squared * par.one_es / (1 - par.es * rh_squared)));
                        return;
                    }

                    if (this->m_proj_parm.mode == equit) {
                        // Equatorial case. Forward case equations can be simplified as:
                        // x = nu * cosphi * sinlam    =>  x^2 * (1 - es * sinphi^2) = (1 - sinphi^2) * sinlam^2
                        // y = nu * sinphi * (1 - es)  =>  y^2 / ((1 - es)^2 + y^2 * es) = sinphi^2

                        if (boost::math::pow<2>(xy_x) + boost::math::pow<2>(xy_y) / par.one_es > 1) {
                            BOOST_THROW_EXCEPTION( projection_exception(error_tolerance_condition) );
                        }

                        const T sinphi_squared = (xy_y == 0.0) ? 0.0 : 1.0 / ( boost::math::pow<2>(par.one_es / xy_y) + par.es );
                        if (sinphi_squared >= T{1} - std::numeric_limits<T>::epsilon()) {
                            lp_lat = std::copysign( half_pi, xy_y );
                            lp_lon = 0.0;
                            return;
                        }

                        lp_lat = std::copysign( asin(sqrt(sinphi_squared)), xy_y );
                        const T sinlam = xy_x * sqrt((1 - par.es * sinphi_squared) / (1.0 - sinphi_squared));
                        lp_lon = asin(std::clamp(sinlam, -1.0, 1.0));
                        return;
                    }

                    T rh, cosc, sinc;

                    if ((sinc = (rh = boost::math::hypot(xy_x, xy_y))) > 1.) {
                        if ((sinc - 1.) > epsilon10) {
                            BOOST_THROW_EXCEPTION( projection_exception(error_tolerance_condition) );
                        }
                        sinc = 1.;
                    }
                    cosc = sqrt(1. - sinc * sinc); /* in this range OK */
                    if (fabs(rh) <= epsilon10) {
                        lp_lat = par.phi0;
                        lp_lon = 0.0;
                    } else {
                        lp_lat = cosc * this->m_proj_parm.sinph0 + xy_y * sinc * this->m_proj_parm.cosph0 /rh;
                        const T xy_y_temp = (cosc - this->m_proj_parm.sinph0 * lp_lat) * rh;
                        const T xy_x_temp = xy_x * sinc * this->m_proj_parm.cosph0;
                        if (fabs(lp_lat) >= 1.)
                            lp_lat = lp_lat < 0. ? -half_pi : half_pi;
                        else
                            lp_lat = asin(lp_lat);
                        lp_lon = (xy_y_temp == 0. && (this->m_proj_parm.mode == obliq || this->m_proj_parm.mode == equit))
                             ? (xy_x_temp == 0. ? 0. : xy_x_temp < 0. ? -half_pi : half_pi)
                                           : atan2(xy_x_temp, xy_y_temp);
                    }

                    const T cosphi0 = this->m_proj_parm.cosph0;
                    const T sinphi0 = this->m_proj_parm.sinph0;
                    const T nu0 = this->m_proj_parm.nu0;
                    for (int i = 0; i < 20; i++) {
                        const double cosphi = cos(lp_lat);
                        const double sinphi = sin(lp_lat);
                        const double coslam = cos(lp_lon);
                        const double sinlam = sin(lp_lon);
                        const double one_minus_es_sinphi2 = 1.0 - par.es * sinphi * sinphi;
                        const double nu = 1.0 / sqrt(one_minus_es_sinphi2);
                        T xy_new_x;
                        T xy_new_y;
                        xy_new_x = nu * cosphi * sinlam;
                        xy_new_y = nu * (sinphi * cosphi0 - cosphi * sinphi0 * coslam) +
                                   par.es * (nu0 * sinphi0 - nu * sinphi) * cosphi0;
                        const double rho = par.one_es * nu / one_minus_es_sinphi2;
                        const double J11 = -rho * sinphi * sinlam;
                        const double J12 = nu * cosphi * coslam;
                        const double J21 =
                            rho * (cosphi * cosphi0 + sinphi * sinphi0 * coslam);
                        const double J22 = nu * sinphi0 * cosphi * sinlam;
                        const double D = J11 * J22 - J12 * J21;
                        const double dx = xy_x - xy_new_x;
                        const double dy = xy_y - xy_new_y;
                        const double dphi = (J22 * dx - J12 * dy) / D;
                        const double dlam = (-J21 * dx + J11 * dy) / D;
                        lp_lat += dphi;
                        if (lp_lat > half_pi) {
                            lp_lat = half_pi - (lp_lat - half_pi);
                            lp_lon = adjlon(lp_lon + M_PI);
                        } else if (lp_lat < -half_pi) {
                            lp_lat = -half_pi + (-half_pi - lp_lat);
                            lp_lon = adjlon(lp_lon + M_PI);
                        }
                        lp_lon += dlam;
                        if (fabs(dphi) < 1e-12 && fabs(dlam) < 1e-12) {
                            return;
                        }
                    }
                    BOOST_THROW_EXCEPTION( projection_exception(error_non_convergent) );
                }

                static inline std::string get_name()
                {
                    return "ortho_ellipsoid";
                }

            };

            // Orthographic
            template <typename Parameters, typename T>
            inline void setup_ortho(Parameters& par, par_ortho<T>& proj_parm)
            {
                if (fabs(fabs(par.phi0) - geometry::math::half_pi<T>()) <= epsilon10)
                    proj_parm.mode = par.phi0 < 0. ? s_pole : n_pole;
                else if (fabs(par.phi0) > epsilon10) {
                    proj_parm.mode = obliq;
                    proj_parm.sinph0 = sin(par.phi0);
                    proj_parm.cosph0 = cos(par.phi0);
                    proj_parm.nu0 = 1.0 / sqrt(1.0 - par.es * proj_parm.sinph0 * proj_parm.sinph0);
                } else
                    proj_parm.mode = equit;
            }

    }} // namespace detail::ortho
    #endif // doxygen

    /*!
        \brief Orthographic projection
        \ingroup projections
        \tparam Geographic latlong point type
        \tparam Cartesian xy point type
        \tparam Parameters parameter type
        \par Projection characteristics
         - Azimuthal
         - Spheroid
        \par Example
        \image html ex_ortho.gif
    */
    template <typename T, typename Parameters>
    struct ortho_spheroid : public detail::ortho::base_ortho_spheroid<T, Parameters>
    {
        template <typename Params>
        inline ortho_spheroid(Params const& , Parameters & par)
        {
            detail::ortho::setup_ortho(par, this->m_proj_parm);
        }
    };

    template <typename T, typename Parameters>
    struct ortho_ellipsoid : public detail::ortho::base_ortho_ellipsoid<T, Parameters>
    {
        template <typename Params>
        inline ortho_ellipsoid(Params const& , Parameters & par)
        {
            detail::ortho::setup_ortho(par, this->m_proj_parm);
        }
    };

    #ifndef DOXYGEN_NO_DETAIL
    namespace detail
    {

        // Static projection
        BOOST_GEOMETRY_PROJECTIONS_DETAIL_STATIC_PROJECTION_FI2(srs::spar::proj_ortho, ortho_spheroid, ortho_ellipsoid)

        // Factory entry(s)
        BOOST_GEOMETRY_PROJECTIONS_DETAIL_FACTORY_ENTRY_FI2(ortho_entry, ortho_spheroid, ortho_ellipsoid)

        BOOST_GEOMETRY_PROJECTIONS_DETAIL_FACTORY_INIT_BEGIN(ortho_init)
        {
            BOOST_GEOMETRY_PROJECTIONS_DETAIL_FACTORY_INIT_ENTRY(ortho, ortho_entry)
        }

    } // namespace detail
    #endif // doxygen

} // namespace projections

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_PROJECTIONS_ORTHO_HPP

