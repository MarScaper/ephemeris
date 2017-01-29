/*
 * Ephemeris.hpp
 */
/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Ephemeris_h
#define Ephemeris_h

#include "Calendar.hpp"
#include "VSOP87.h"

/*! This structure describes equatorial coordinates. */
struct EquatorialCoordinates
{
    /*! Floating value for Right Ascension. */
    float ra;
    
    /*! Floating value for Declination */
    float dec;
};

/*! This structure describes horizontal coordinates. */
struct HorizontalCoordinates
{
    /*! Floating value for altitude. */
    float alt;
    
    /*! Floating value for azimuth */
    float azi;
};

/*! This structure describes Heliocentric ecliptic coordinates. */
struct HeliocentricCoordinates
{
    /*! Floating value for ecliptic longitude. */
    float lon;
    
    /*! Floating value for ecliptic latitude.*/
    float lat;
    
    /*! Floating value for radius vector (distance from Sun). */
    float radius;
};

/*! This structure describes geocentric coordinates. */
struct GeocentricCoordinates
{
    /*! Floating value for longitude. */
    float lon;
    
    /*! Floating value for latitude.*/
    float lat;
};

/*! This structure describes rectangular coordinates. */
struct RectangularCoordinates
{
    float x;
    float y;
    float z;
};

/*! This structure describes available solar system objects for computation of ephemerides. */
enum SolarSystemObjectIndex
{
    Sun        = 0,
    Mercury    = 1,
    Venus      = 2,
    Earth      = 3,
    Mars       = 4,
    Jupiter    = 5,
    Saturn     = 6,
    Uranus     = 7,
    Neptune    = 8,
    
    // Not implemented yet 
    //EarthsMoon = 9
};

/*! This structure describes a planet for a specific date and time. */
struct SolarSystemObject
{
    /*! Equatorial coordinates (RA/Dec). */
    EquatorialCoordinates   equaCoordinates;
    
    /*! Horizontal coordinates (Alt/Az). */
    HorizontalCoordinates horiCoordinates;
    
    /*! Apparent diameter from earth in arc seconds. */
    float diameter;
    
    /*! Distance from earth in astronomical unit. */
    float distance;
};

/*! This structure describes planetary orbit. */
struct PlanetayOrbit
{
    /*! Mean longitude. */
    float L;
    
    /*! Semimajor axis. */
    float a;
    
    /*! Eccentricity. */
    float e;
    
    /*! Inclination. */
    float i;
    
    /*! Longitude ascending node. */
    float omega;
    
    /*! Perihelion. */
    float pi;
    
    /*! Mean anomaly. */
    float M;
    
    /*! Perihelion argument. */
    float w;
};

/*!
 * This class is used for astronomical calculations. The code is based on the book "Astronomical Algorithms" by Jean Meeus.
 */
class Ephemeris
{
    
public:
    
    /*! Set location on earth (used for horizontal coordinates conversion). */
    static void setLocationOnEarth(float floatingLatitude, float floatingLongitude);
    
    /*! Set location on earth (used for horizontal coordinates conversion). */
    static void setLocationOnEarth(float latDegrees, float latMinutes, float latSeconds,
                                   float lonDegrees, float lonMinutes, float lonSeconds);
    
    
    
    /*! Convert floating hours to integer hours, minutes, seconds. */
    static void  floatingHoursToHoursMinutesSeconds(float floatingHours, int *hours, int *minutes, float *seconds);
    
    /*! Convert integer hours, minutes, seconds to floating hours. */
    static float hoursMinutesSecondsToFloatingHours(int hours, int minutes, float seconds);
    
    /*! Convert floating degrees to integer degrees, minutes, seconds. */
    static void  floatingDegreesToDegreesMinutesSeconds(float floatingDegrees, int *degrees, int *minutes, float *seconds);
    
    /*! Convert integer degrees, minutes, seconds to floating degrees. */
    static float degreesMinutesSecondsToFloatingDegrees(int degrees, int minutes, float seconds);
    
    
    
    /*! Compute sun coordinates in the sky (R.A.,Dec) for a specific date and time.
     *  Reference: Chapter 16, page 63: Les coordonnées du soleil. */
    static EquatorialCoordinates sunForDateAndTime(unsigned int day,  unsigned int month,  unsigned int year,
                                                   unsigned int hours, unsigned int minutes, unsigned int seconds);
    
    /*! Compute solar system object for a specific date, time and location on earth (if location has been initialized first). */
    static SolarSystemObject solarSystemObjectAtDateAndTime(SolarSystemObjectIndex planet,
                                                            unsigned int day,  unsigned int month,  unsigned int year,
                                                            unsigned int hours, unsigned int minutes, unsigned int seconds);
    
private:
    
    /*! Compute apparent sideral time (in floating hours) for a given date and time.
     *  Reference: Chapter 7, page 35: Temps sidéral à Greenwich. */
    static float apparentSideralTime(unsigned int day,   unsigned int month,   unsigned int year,
                                     unsigned int hours, unsigned int minutes, unsigned int seconds);
    
    /*! Compute heliocentric coordinates. */
    static HeliocentricCoordinates heliocentricCoordinatesForPlanetAndT(SolarSystemObjectIndex planet, float T);
    
    /*! Compute Kepler equation. */
    static float kepler(float M, float e);
    
    /*! Convert equatorial coordinates to horizontal coordinates. */
    static HorizontalCoordinates equatorialToHorizontal(float H, float delta, float phi);
    
    /*! Convert ecliptic coordinates to equatorial coordinates. */
    static EquatorialCoordinates EclipticToEquatorial(float lambda, float beta, float epsilon);
    
    /*! Convert heliocentric coordinates to rectangular coordinates. */
    static RectangularCoordinates HeliocentricToRectangular(HeliocentricCoordinates hc, HeliocentricCoordinates hc0);
    
    static float meanGreenwichSiderealTimeAtDateAndTime(unsigned int day,   unsigned int month,   unsigned int year,
                                                        unsigned int hours, unsigned int minutes, unsigned int seconds);
    
    /*! Compute the true obliquity (angle in floating degrees) of the ecliptic,
     *  obliquity, delta obliquity and delta nutation for T.
     *  Reference: Chapter 13, page 53: Nutation et obliquité de l'écliptique. */
    static float obliquityAndNutationForT(float T, float *deltaObliquity, float *deltaNutation);
    
    /*! Compute planet informations for T. */
    static PlanetayOrbit planetayOrbitForPlanetAndT(SolarSystemObjectIndex planet, float T);
    
    static EquatorialCoordinates equatorialCoordinatesForSunAtJD(JulianDay jd, float *distance, GeocentricCoordinates *gCoordinates);
    
    /*! Compute equatorial coordinates (and geocentric if needed) for a a specific Julian day. */
    static EquatorialCoordinates equatorialCoordinatesForPlanetAtJD(SolarSystemObjectIndex planet, JulianDay jd,
                                                                    float *distance, GeocentricCoordinates *gCoordinates);
    
    /*! Compute VSOP87 coefficients for T. */
    static float sumVSOP87Coefs(const VSOP87Coefficient *valuePlanetCoefficients, int coefCount, float T);
    
};

#endif
