/*
 * Ephemeris.cpp
 *
 * Copyright (c) 2017 by Sebastien MARCHAND (Web:www.marscaper.com, Email:sebastien@marscaper.com)
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

#if ARDUINO
#include <Arduino.h>
#endif
#include <stdio.h>
#include <math.h>

#include "Ephemeris.hpp"


#ifndef PI
#define PI 3.1415926535
#endif

// Trigonometry using degrees
#define SIND(value)   sin((value)*0.0174532925)
#define COSD(value)   cos((value)*0.0174532925)
#define TAND(value)   tan((value)*0.0174532925)

#define ACOSD(value)  (acos((value))*57.2957795131);
#define ATAND(value)  (atan((value))*57.2957795131);

// Limit range
#define LIMIT_DEGREES_TO_360(value) (value) >= 0 ? ((value)-(long)((value)*0.0027777778)*360) : (((value)-(long)((value)*0.0027777778)*360)+360)
#define LIMIT_HOURS_TO_24(value) (value) >= 0 ? ((value)-(long)((value)*0.0416666667)*24) : ((value)+24)
#define LIMIT_DEC_TO_90(value) (value) 

// Convert degrees
#define DEGREES_TO_RADIANS(value) ((value)*0.0174532925)
#define DEGREES_TO_HOURS(value) ((value)*0.0666666667)
#define DEGREES_MINUTES_SECONDES_TO_SECONDS(deg,min,sec) ((FLOAT)(deg)*3600+(FLOAT)(min)*60+(FLOAT)sec)
#define DEGREES_MINUTES_SECONDS_TO_DECIMAL_DEGREES(deg,min,sec) (deg) >= 0 ? ((FLOAT)(deg)+(FLOAT)(min)*0.0166666667+(FLOAT)(sec)*0.0002777778) : ((FLOAT)(deg)-(FLOAT)(min)*0.0166666667-(FLOAT)(sec)*0.0002777778)

// Convert radians
#define RADIANS_TO_DEGREES(value) ((value)*57.2957795131)
#define RADIANS_TO_HOURS(value) ((value)*3.81971863)

// Convert hours
#define HOURS_TO_RADIANS(value) ((value)*0.261799388)
#define HOURS_MINUTES_SECONDS_TO_SECONDS(hour,min,sec) ((FLOAT)(hour)*3600+(FLOAT)(min)*60+(FLOAT)sec)
#define HOURS_MINUTES_SECONDS_TO_DECIMAL_HOURS(hour,min,sec) ((FLOAT)(hour)+(FLOAT)(min)*0.0166666667+(FLOAT)(sec)*0.0002777778)
#define HOURS_TO_DEGREES(value) ((value)*15)

// Convert seconds
#define SECONDS_TO_DECIMAL_DEGREES(value) ((FLOAT)(value)*0.0002777778)
#define SECONDS_TO_DECIMAL_HOURS(value) ((FLOAT)(value)*0.0002777778)

#define T_WITH_JD(day,time)((day-2451545.0+time)/36525)

// Observer's coordinates on Earth
static FLOAT latitudeOnEarth      = NAN;
static FLOAT longitudeOnEarth     = NAN;
static FLOAT longitudeOnEarthSign = -1;
static int   altitudeOnEarth      = NAN;

void Ephemeris::floatingHoursToHoursMinutesSeconds(FLOAT floatingHours, int *hours, int *minutes, FLOAT *seconds)
{
    floatingHours = LIMIT_HOURS_TO_24(floatingHours);
    
    // Calculate hours, minutes, seconds
    *hours   = (int)floatingHours;
    *minutes = floatingHours*60-(long)(*hours)*60;
    *seconds = floatingHours*3600-(long)(*hours)*3600-(long)(*minutes)*60;
}

FLOAT Ephemeris::hoursMinutesSecondsToFloatingHours(int degrees, int minutes, FLOAT seconds)
{
    return HOURS_MINUTES_SECONDS_TO_DECIMAL_HOURS(degrees, minutes, seconds);
}

void Ephemeris::floatingDegreesToDegreesMinutesSeconds(FLOAT floatintDegrees, int *degree, int *minute, FLOAT *second)
{
    // Calculate hour,minute,second
    if( floatintDegrees >= 0 )
    {
        *degree   = (unsigned int)floatintDegrees;
        *minute = floatintDegrees*60-(long)(*degree*60);
        *second = floatintDegrees*3600-(long)(*degree*3600)-(long)(*minute*60);
    }
    else
    {
        floatintDegrees = floatintDegrees*-1;
        
        *degree   = (unsigned int)floatintDegrees;
        *minute = floatintDegrees*60-(long)(*degree*60);
        *second = floatintDegrees*3600-(long)(*degree*3600)-(long)(*minute*60);
        
        *degree = *degree*-1;
    }
    
    return;
}

FLOAT Ephemeris::degreesMinutesSecondsToFloatingDegrees(int degrees, int minutes, FLOAT seconds)
{
    return DEGREES_MINUTES_SECONDS_TO_DECIMAL_DEGREES(degrees,minutes,seconds);
}

HorizontalCoordinates Ephemeris::equatorialToHorizontalCoordinatesAtDateAndTime(EquatorialCoordinates eqCoordinates,
                                                                                unsigned int day,  unsigned int month,  unsigned int year,
                                                                                unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    HorizontalCoordinates hCoordinates;
    
    if( !isnan(longitudeOnEarth) && !isnan(latitudeOnEarth) )
    {
        JulianDay jd = Calendar::julianDayForDate(day, month, year);
        
        FLOAT T = T_WITH_JD(jd.day,jd.time);
        
        FLOAT meanSideralTime = meanGreenwichSiderealTimeAtDateAndTime(day, month, year, hours, minutes, seconds);
        
        FLOAT deltaNutation;
        FLOAT epsilon = obliquityAndNutationForT(T, NULL, &deltaNutation);
        
        // Apparent sideral time in floating hours
        FLOAT theta0 = meanSideralTime + (deltaNutation/15*COSD(epsilon))/3600;
        
        
        // Geographic longitude in floating hours
        FLOAT L = DEGREES_TO_HOURS(longitudeOnEarth*longitudeOnEarthSign);
        
        // Geographic latitude in floating degrees
        FLOAT phi = latitudeOnEarth;
        
        // Local angle in floating degrees
        FLOAT H = (theta0-L-eqCoordinates.ra)*15;
        
        hCoordinates = equatorialToHorizontal(H,eqCoordinates.dec,phi);
    }
    else
    {
        hCoordinates.alt = NAN;
        hCoordinates.azi = NAN;
    }
    
    return hCoordinates;
}

EquatorialCoordinates Ephemeris::horizontalToEquatorialCoordinatesAtDateAndTime(HorizontalCoordinates hCoordinates,
                                                                                unsigned int day,  unsigned int month,  unsigned int year,
                                                                                unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    EquatorialCoordinates eqCoordinates;
    
    if( !isnan(longitudeOnEarth) && !isnan(latitudeOnEarth) )
    {
        
        JulianDay jd = Calendar::julianDayForDate(day, month, year);
        
        FLOAT T = T_WITH_JD(jd.day,jd.time);
        
        FLOAT meanSideralTime = meanGreenwichSiderealTimeAtDateAndTime(day, month, year, hours, minutes, seconds);
        
        FLOAT deltaNutation;
        FLOAT epsilon = obliquityAndNutationForT(T, NULL, &deltaNutation);
        
        // Apparent sideral time in floating hours
        FLOAT theta0 = meanSideralTime + (deltaNutation/15*COSD(epsilon))/3600;
        
        
        // Geographic longitude in floating hours
        FLOAT L = DEGREES_TO_HOURS(longitudeOnEarth*longitudeOnEarthSign);
        
        // Geographic latitude in floating degrees
        FLOAT phi = latitudeOnEarth;
        
        // Compute local horizontal coordinates (note: RA will contain H)
        eqCoordinates = horizontalToEquatorial(hCoordinates.azi, hCoordinates.alt, phi);
        
        // Compute RA according to H
        // RA = theta0 - L - H;
        eqCoordinates.ra = theta0 - L - eqCoordinates.ra;
        eqCoordinates.ra = LIMIT_HOURS_TO_24(eqCoordinates.ra);
    }
    else
    {
        eqCoordinates.ra  = NAN;
        eqCoordinates.dec = NAN;
    }
    
    return eqCoordinates;
}

FLOAT Ephemeris::apparentSideralTime(unsigned int day,  unsigned int month,  unsigned int year,
                                     unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    JulianDay jd = Calendar::julianDayForDate(day, month, year);
    
    FLOAT T        = T_WITH_JD(jd.day,jd.time);
    FLOAT TSquared = T*T;
    FLOAT TCubed   = TSquared*T;
    
    FLOAT theta0 = 100.46061837 + T*36000.770053608 + TSquared*0.000387933  - TCubed/38710000;
    theta0 = LIMIT_DEGREES_TO_360(theta0);
    theta0 = DEGREES_TO_HOURS(theta0);
    
    FLOAT time = HOURS_MINUTES_SECONDS_TO_DECIMAL_HOURS(hours,minutes,seconds);
    
    FLOAT apparentSideralTime = theta0 + 1.00273790935 * time;
    apparentSideralTime = LIMIT_HOURS_TO_24(apparentSideralTime);
    
    return apparentSideralTime;
}

FLOAT Ephemeris::obliquityAndNutationForT(FLOAT T, FLOAT *deltaObliquity, FLOAT *deltaNutation)
{
    FLOAT TSquared = T*T;
    FLOAT TCubed   = TSquared*T;
    
    FLOAT Ls = 280.4565 + T*36000.7698   + TSquared*0.000303;
    Ls = LIMIT_DEGREES_TO_360(Ls);
    
    FLOAT Lm = 218.3164 + T*481267.8812  - TSquared*0.001599;
    Lm = LIMIT_DEGREES_TO_360(Lm);
    
    FLOAT Ms = 357.5291 + T*35999.0503   - TSquared*0.000154;
    Ms = LIMIT_DEGREES_TO_360(Ms);
    
    FLOAT Mm = 134.9634 + T*477198.8675  + TSquared*0.008721;
    Mm = LIMIT_DEGREES_TO_360(Mm);
    
    FLOAT omega = 125.0443 - T*1934.1363 + TSquared*0.008721;
    omega = LIMIT_DEGREES_TO_360(omega);
    
    // Delta Phi
    FLOAT dNutation =
    -(17.1996 + 0.01742*T) * SIND(omega)
    -(1.3187  + 0.00016*T) * SIND(2*Ls)
    - 0.2274               * SIND(2*Lm)
    + 0.2062               * SIND(2*omega)
    +(0.1426  - 0.00034*T) * SIND(Ms)
    + 0.0712               * SIND(Mm)
    -(0.0517  - 0.00012*T) * SIND(2*Ls+Ms)
    - 0.0386               * SIND(2*Lm-omega)
    - 0.0301               * SIND(2*Lm+Mm)
    + 0.0217               * SIND(2*Ls-Ms)
    - 0.0158               * SIND(2*Ls-2*Lm+Mm)
    + 0.0129               * SIND(2*Ls-omega)
    + 0.0123               * SIND(2*Lm-Mm);
    
    if( deltaNutation )
    {
        *deltaNutation = dNutation;
    }
    
    // Delta Eps
    FLOAT dObliquity =
    +(9.2025  + 0.00089*T) * COSD(omega)
    +(0.5736  - 0.00031*T) * COSD(2*Ls)
    + 0.0977               * COSD(2*Lm)
    - 0.0895               * COSD(2*omega)
    + 0.0224               * COSD(2*Ls+Ms)
    + 0.0200               * COSD(2*Lm-omega)
    + 0.0129               * COSD(2*Lm + Mm)
    - 0.0095               * COSD(2*Ls-Ms)
    - 0.0070               * COSD(2*Ls-omega);
    
    if( deltaObliquity )
    {
        *deltaObliquity = dObliquity;
    }
    
    FLOAT eps0 = DEGREES_MINUTES_SECONDES_TO_SECONDS(23,26,21.448)-T*46.8150-TSquared*0.00059+TCubed*0.001813;
    
    FLOAT obliquity = eps0 + dObliquity;
    obliquity = SECONDS_TO_DECIMAL_DEGREES(obliquity);
    
    return obliquity;
}

#if !DISABLE_PLANETS
FLOAT Ephemeris::sumELP2000Coefs(const FLOAT *moonMultCoefficients, const ELP2000Coefficient *moonAngleCoefficients, int coefCount,
                                 FLOAT E, FLOAT D, FLOAT M, FLOAT Mp, FLOAT F, bool squareMultiplicator)
{
    // Parse each value in coef table
    FLOAT value = 0;
    for(int numCoef=0; numCoef<coefCount; numCoef++)
    {
        // Get coef
        FLOAT multiplicator;
        ELP2000Coefficient moonAngle;;
        
#if ARDUINO
        // We limit SRAM usage by using flash memory (PROGMEM)
        memcpy_P(&multiplicator, &moonMultCoefficients[numCoef],  sizeof(FLOAT));
        memcpy_P(&moonAngle,     &moonAngleCoefficients[numCoef], sizeof(ELP2000Coefficient));
#else
        multiplicator = moonMultCoefficients[numCoef];
        moonAngle     = moonAngleCoefficients[numCoef];
#endif
        
        FLOAT angle =  moonAngle.D*D + moonAngle.M*M + moonAngle.Mp*Mp + moonAngle.F*F;
        
        FLOAT res;
        if( moonMultCoefficients != RMoonCoefficients )
        {
            res  = multiplicator * SIND(angle);
        }
        else
        {
            // R coef use cos and not sin.
            res  = multiplicator * COSD(angle);
        }
        
        if(squareMultiplicator)
        {
            // To avoid out of range issue with single precision
            // we've stored sqrt(value) and not the value. As a result we need to square it back.
            res *= multiplicator;
            
            // Now respect original sign
            if( multiplicator < 0 )
            {
                res *= -1;
            }
        }
        
        if( moonAngle.Ec == 1 )
        {
            // E
            res *= E;
        }
        else if( moonAngle.Ec == 2 )
        {
            // E^2
            res *= E;
            res *= E;
        }
        
        value += res;
    }
    
    return value;
}
#endif

#if !DISABLE_PLANETS
EquatorialCoordinates  Ephemeris::equatorialCoordinatesForEarthsMoonAtJD(JulianDay jd, FLOAT *distance)
{
    FLOAT T        = T_WITH_JD(jd.day,jd.time);
    FLOAT TSquared = T*T;
    FLOAT TCubed   = TSquared*T;
    FLOAT T4       = TCubed*T;
    
    // Mean Longitude of the Moon
    FLOAT Lp = 218.31644735 + 481267.88122838*T - 0.001579944*TSquared + TCubed/538841 - T4/538841;
    Lp = LIMIT_DEGREES_TO_360(Lp);
    
    // Mean anomaly of the Sun
    FLOAT M = 357.5291092 + T*35999.0502909  - TSquared*0.0001536;
    M = LIMIT_DEGREES_TO_360(M);
    
    // Mean anomaly of Moon
    FLOAT Mp = 134.96339622 + 477198.86750067*T + 0.00872053*TSquared + TCubed/69699;
    Mp = LIMIT_DEGREES_TO_360(Mp);
    
    // Mean elongation of Moon
    FLOAT D = 297.85019172 + 445267.11139756*T - 0.00190272*TSquared + TCubed/545868;
    D = LIMIT_DEGREES_TO_360(D);
    
    // Mean distance at ascending node
    FLOAT F = 93.27209769 + 483202.01756053*T - 0.00367481*TSquared - TCubed/3525955;
    F = LIMIT_DEGREES_TO_360(F);
    
    // A1, A2, A3
    FLOAT A1 = 119.75 + 131.849*T;
    A1 = LIMIT_DEGREES_TO_360(A1);
    FLOAT A2 =  53.09 + 479264.290*T;
    A2 = LIMIT_DEGREES_TO_360(A2);
    FLOAT A3 = 313.45 + 481266.484*T;
    A3 = LIMIT_DEGREES_TO_360(A3);
    
    
    FLOAT E = 1 - 0.002516*T - 0.0000074*TCubed;
    E = LIMIT_DEGREES_TO_360(E);
    
    
    FLOAT sumL = sumELP2000Coefs( LMoonCoefficients, LMoonAngleCoefficients, sizeof(LMoonCoefficients)/sizeof(FLOAT), E, D, M, Mp, F, false);
    sumL += 3958*SIND(A1) + 1962*SIND(Lp-F) + 318*SIND(A2);
    
    FLOAT sumB = sumELP2000Coefs(BMoonCoefficients, BMoonAngleCoefficients, sizeof(BMoonCoefficients)/sizeof(FLOAT), E, D, M, Mp, F, false);
    sumB += -2235*SIND(Lp) + 382*SIND(A3) + 175*SIND(A1-F) + 175*SIND(A1+F) + 127*SIND(Lp-Mp) - 115*sin(Lp+Mp);
    
    FLOAT sumR = sumELP2000Coefs(RMoonCoefficients, RMoonAngleCoefficients, sizeof(RMoonCoefficients)/sizeof(FLOAT), E, D, M, Mp, F, true);
    
    // Geocentic longitude
    FLOAT lambda = Lp + sumL/1000000;   // Degrees
    
    // Geocentric latitude
    FLOAT beta   = sumB/1000000;        // Degrees
    
    // Distance
    FLOAT delta  = 385000.56+sumR/1000; // Kilometers
    
    // Convert kilometers to AU
    FLOAT dist = delta*6.68459e-9;
    
    if( distance )
    {
        *distance = dist;
    }
    
    // Obliquity and Nutation
    FLOAT deltaNutation;
    FLOAT epsilon = Ephemeris::obliquityAndNutationForT(T, NULL, &deltaNutation);
    
    // Intergrate nutation
    lambda += deltaNutation/3600;
    
    EquatorialCoordinates eqCoord = EclipticToEquatorial(lambda, beta, epsilon);
    
    //
    // Geocentric to topocentric conversion
    //
    
    FLOAT paralax = 8.794/(dist)/3600.f;
    
    FLOAT meanSideralTime = meanGreenwichSiderealTimeAtJD(jd);
    
    // Geographic longitude in floating hours
    FLOAT L = DEGREES_TO_HOURS(longitudeOnEarth*longitudeOnEarthSign);
    
    // Apparent sideral time in floating hours
    FLOAT theta0 = meanSideralTime + (deltaNutation/15*COSD(epsilon))/3600;
    
    // Local angle in floating degrees
    FLOAT H = (theta0-L-eqCoord.ra)*15;
    
    // Geocentric rectangular coordinates of the observer (Chap 6 p33)
    FLOAT u = ATAND(0.99664719*TAND(latitudeOnEarth));
    FLOAT rhoSinDeltaPrime    = 0.99664719*SIND(u)+(altitudeOnEarth/6378140)*SIND(latitudeOnEarth);
    FLOAT rhoCosDeltaPrime = COSD(u)+(altitudeOnEarth/6378140)*COSD(latitudeOnEarth);
    
    // Paralax correction (Chap 27 p 105)
    FLOAT detlaRa = ATAND(-rhoCosDeltaPrime*SIND(paralax)*SIND(H)/(COSD(eqCoord.dec)-rhoCosDeltaPrime*SIND(paralax)*COSD(H)));
    FLOAT detlaDec = ATAND((SIND(eqCoord.dec)-rhoSinDeltaPrime*SIND(paralax))*COSD(detlaRa)/(COSD(eqCoord.dec)-rhoCosDeltaPrime*SIND(paralax)*COSD(H)));
    
    eqCoord.ra += detlaRa/15.f;
    eqCoord.dec = detlaDec;

    return eqCoord;
}
#endif

#if !DISABLE_PLANETS
EquatorialCoordinates  Ephemeris::equatorialCoordinatesForSunAtJD(JulianDay jd, FLOAT *distance)
{
    EquatorialCoordinates sunCoordinates;
    
    FLOAT T        = T_WITH_JD(jd.day,jd.time);
    FLOAT TSquared = T*T;
    
    FLOAT L0 = 280.46646 + T*36000.76983 + TSquared*0.0003032;
    L0 = LIMIT_DEGREES_TO_360(L0);
    
    FLOAT M = 357.5291092 + T*35999.0502909  - TSquared*0.0001536;
    M = LIMIT_DEGREES_TO_360(M);
    
    FLOAT e = 0.016708634 - T*0.000042037 - TSquared*0.0000001267;
    
    FLOAT C =
    +(1.914602 - T*0.004817 - TSquared*0.000014) * SIND(M)
    +(0.019993 - T*0.000101                    ) * SIND(2*M)
    + 0.000289                                   * SIND(3*M);
    
    FLOAT O = L0 + C;
    
    FLOAT v = M  + C;
    
    // Improved precision for O according to page 65
    {
        FLOAT Av = 351.52 + 22518.4428*T;  // Mars
        FLOAT Bv = 253.14 + 45036.8857*T;  // Venus
        FLOAT Cj = 157.23 + 32964.4673*T;  // Jupiter
        FLOAT Dm = 297.85 + 445267.1117*T; // Moon
        FLOAT E  = 252.08 + 20.19 *T;
        FLOAT H  = 42.43  + 65928.9358*T;
        
        O +=
        + 0.00134 * COSD(Av)
        + 0.00153 * COSD(Bv)
        + 0.00200 * COSD(Cj)
        + 0.00180 * SIND(Dm)
        + 0.00196 * SIND(E);
        
        v +=
        + 0.00000542 * SIND(Av)
        + 0.00001576 * SIND(Bv)
        + 0.00001628 * SIND(Cj)
        + 0.00003084 * COSD(Dm)
        + 0.00000925 * SIND(H);
    }
    
    // R
    FLOAT dist = (1.000001018*(1-e*e))/(1+e*COSD(v));
    if( distance ) *distance = dist;

    FLOAT omega = 125.04 - 1934.136*T;
    
    FLOAT lambda = O - 0.00569 - 0.00478 * SIND(omega);
    
    FLOAT deltaNutation;
    FLOAT epsilon = obliquityAndNutationForT(T, NULL, &deltaNutation);
    
    FLOAT eps = epsilon + 0.00256*COSD(omega);
    
    // Alpha   (Hour=Deg/15.0)
    sunCoordinates.ra = atan2(SIND(lambda)*COSD(eps),COSD(lambda));
    sunCoordinates.ra = RADIANS_TO_HOURS(sunCoordinates.ra);
    sunCoordinates.ra = LIMIT_HOURS_TO_24(sunCoordinates.ra);
    
    // Delta
    sunCoordinates.dec = asin(SIND(eps)*SIND(lambda))*180.0/PI;
    
    EquatorialCoordinates eqCoord = sunCoordinates;
    
    //
    // Geocentric to topocentric conversion
    //
    
    FLOAT paralax = 8.794/(dist)/3600.f;
    
    FLOAT meanSideralTime = meanGreenwichSiderealTimeAtJD(jd);
    
    // Geographic longitude in floating hours
    FLOAT L = DEGREES_TO_HOURS(longitudeOnEarth*longitudeOnEarthSign);
    
    // Apparent sideral time in floating hours
    FLOAT theta0 = meanSideralTime + (deltaNutation/15*COSD(epsilon))/3600;
    
    // Local angle in floating degrees
    FLOAT H = (theta0-L-eqCoord.ra)*15;
    
    // Geocentric rectangular coordinates of the observer (Chap 6 p33)
    FLOAT u = ATAND(0.99664719*TAND(latitudeOnEarth));
    FLOAT rhoSinDeltaPrime    = 0.99664719*SIND(u)+(altitudeOnEarth/6378140)*SIND(latitudeOnEarth);
    FLOAT rhoCosDeltaPrime = COSD(u)+(altitudeOnEarth/6378140)*COSD(latitudeOnEarth);
    
    // Paralax correction (Chap 27 p 105)
    FLOAT detlaRa = ATAND(-rhoCosDeltaPrime*SIND(paralax)*SIND(H)/(COSD(eqCoord.dec)-rhoCosDeltaPrime*SIND(paralax)*COSD(H)));
    FLOAT detlaDec = ATAND((SIND(eqCoord.dec)-rhoSinDeltaPrime*SIND(paralax))*COSD(detlaRa)/(COSD(eqCoord.dec)-rhoCosDeltaPrime*SIND(paralax)*COSD(H)));
    
    eqCoord.ra += detlaRa/15.f;
    eqCoord.dec = detlaDec;
    
    return eqCoord;
}
#endif

#if !DISABLE_PLANETS
PlanetayOrbit Ephemeris::planetayOrbitForPlanetAndT(SolarSystemObjectIndex solarSystemObjectIndex, FLOAT T)
{
    PlanetayOrbit planetayOrbit;
    
    FLOAT TSquared = T*T;
    FLOAT TCubed   = TSquared*T;
    
    switch (solarSystemObjectIndex)
    {
        case Mercury:
            planetayOrbit.L     = 252.250906   + 149474.0722491*T + 0.00030350*TSquared   + 0.000000018*TCubed;
            planetayOrbit.a     = 0.387098310;
            planetayOrbit.e     = 0.20563175   + 0.000020407*T    - 0.0000000283*TSquared - 0.00000000018*TCubed;
            planetayOrbit.i     = 7.004986     + 0.0018215*T      - 0.00001810*TSquared   + 0.000000056*TCubed;
            planetayOrbit.omega = 48.330893    + 1.1861883*T      + 0.00017542*TSquared   + 0.000000215*TCubed;
            planetayOrbit.pi    = 77.456119    + 1.5564776*T      + 0.00029544*TSquared   + 0.000000009*TCubed;
            break;
            
        case Venus:
            planetayOrbit.L     = 181.979801   + 58519.2130302*T + 0.00031014*TSquared   + 0.000000015*TCubed;
            planetayOrbit.a     = 0.723329820;
            planetayOrbit.e     = 0.00677192   - 0.000047765*T   + 0.0000000981*TSquared + 0.00000000046*TCubed;
            planetayOrbit.i     = 3.394662     + 0.0010037*T     - 0.00000088*TSquared   - 0.000000007*TCubed;
            planetayOrbit.omega = 76.679920    + 0.9011206*T     + 0.00040618*TSquared   - 0.000000093*TCubed;
            planetayOrbit.pi    = 131.563703   + 1.4022288*T     - 0.00107618*TSquared   - 0.000005678*TCubed;
            break;
            
        case Earth:
            planetayOrbit.L     = 100.466457   + 36000.7698278*T + 0.00030322*TSquared   + 0.000000020*TCubed;
            planetayOrbit.a     = 1.000001018;
            planetayOrbit.e     = 0.01670863   - 0.000042037*T   - 0.0000001267*TSquared + 0.00000000014*TCubed;
            planetayOrbit.i     = 0;
            planetayOrbit.omega = NAN;
            planetayOrbit.pi    = 102.937348   + 1.17195366*T    + 0.00045688*TSquared   - 0.000000018*TCubed;
            break;
            
        case Mars:
            planetayOrbit.L     = 355.433000   + 19141.6964471*T + 0.00031052*TSquared   + 0.000000016*TCubed;
            planetayOrbit.a     = 1.523679342;
            planetayOrbit.e     = 0.09340065   + 0.000090484*T   - 0.0000000806*TSquared - 0.00000000025*TCubed;
            planetayOrbit.i     = 1.849726     - 0.0006011*T     + 0.00001276*TSquared   - 0.000000007*TCubed;
            planetayOrbit.omega = 49.588093    + 0.7720959*T     + 0.00001557*TSquared   + 0.000002267*TCubed;
            planetayOrbit.pi    = 336.060234   + 1.8410449*T     + 0.00013477*TSquared   + 0.000000536*TCubed;
            break;
            
        case Jupiter:
            planetayOrbit.L     = 34.351519   + 3036.3027748*T  + 0.00022330*TSquared   + 0.000000037*TCubed;
            planetayOrbit.a     = 5.202603209 + 0.0000001913*T;
            planetayOrbit.e     = 0.04849793  + 0.000163225*T   - 0.0000004714*TSquared - 0.00000000201*TCubed;
            planetayOrbit.i     = 1.303267    - 0.0054965*T     + 0.00000466*TSquared   - 0.000000002*TCubed;
            planetayOrbit.omega = 100.464407  + 1.0209774*T     + 0.00040315*TSquared   + 0.000000404*TCubed;
            planetayOrbit.pi    = 14.331207   + 1.6126352*T     + 0.00103042*TSquared   - 0.000004464*TCubed;
            break;
            
        case Saturn:
            planetayOrbit.L     = 50.077444   + 1223.5110686*T + 0.00051908*TSquared   - 0.000000030*TCubed;
            planetayOrbit.a     = 9.554909192 - 0.0000021390*T + 0.000000004*TSquared;
            planetayOrbit.e     = 0.05554814  - 0.0003446641*T - 0.0000006436*TSquared + 0.00000000340*TCubed;
            planetayOrbit.i     = 2.488879    - 0.0037362*T    - 0.00001519*TSquared   + 0.000000087*TCubed;
            planetayOrbit.omega = 113.665503  + 0.8770880*T    - 0.00012176*TSquared   - 0.000002249*TCubed;
            planetayOrbit.pi    = 93.057237   + 1.9637613*T    + 0.00083753*TSquared   + 0.000004928*TCubed;
            break;
            
        case Uranus:
            planetayOrbit.L     = 314.055005   + 429.8640561*T  + 0.00030390*TSquared     + 0.000000026*TCubed;
            planetayOrbit.a     = 19.218446062 - 0.0000000372*T + 0.00000000098*TSquared;
            planetayOrbit.e     = 0.04638122   - 0.000027293*T  + 0.0000000789*TSquared   + 0.00000000024*TCubed;
            planetayOrbit.i     = 0.773197     + 0.0007744*T    + 0.00003749*TSquared     - 0.000000092*TCubed;
            planetayOrbit.omega = 74.005957    + 0.5211278*T    + 0.00133947*TSquared     + 0.000018484*TCubed;
            planetayOrbit.pi    = 173.005291   + 1.4863790*T    + 0.00021406*TSquared     + 0.000000434*TCubed;
            break;
            
        case Neptune:
            planetayOrbit.L     = 304.348665   + 219.8833092*T  + 0.00030882*TSquared     + 0.000000018*TCubed;
            planetayOrbit.a     = 30.110386869 - 0.0000001663*T + 0.00000000069*TSquared;
            planetayOrbit.e     = 0.00945575   + 0.000006033*T                            - 0.00000000005*TCubed;
            planetayOrbit.i     = 1.769953     - 0.0093082*T    - 0.00000708*TSquared     + 0.000000027*TCubed;
            planetayOrbit.omega = 131.784057   + 1.1022039*T    + 0.00025952*TSquared     - 0.000000637*TCubed;
            planetayOrbit.pi    = 48.120276    + 1.4262957*T    + 0.00038434*TSquared     + 0.000000020*TCubed;
            break;
            
        default:
            // Unknow planet
            break;
    }
    
    // Apply limitations
    planetayOrbit.L     = LIMIT_DEGREES_TO_360(planetayOrbit.L);
    planetayOrbit.i     = LIMIT_DEGREES_TO_360(planetayOrbit.i);
    planetayOrbit.omega = LIMIT_DEGREES_TO_360(planetayOrbit.omega);
    planetayOrbit.pi    = LIMIT_DEGREES_TO_360(planetayOrbit.pi);
    
    // Mean anomaly
    planetayOrbit.M = planetayOrbit.L - planetayOrbit.pi;
    planetayOrbit.M = LIMIT_DEGREES_TO_360(planetayOrbit.M);
    
    planetayOrbit.w = planetayOrbit.pi-planetayOrbit.omega;
    planetayOrbit.w = LIMIT_DEGREES_TO_360(planetayOrbit.w);
    
    return planetayOrbit;
}
#endif

FLOAT Ephemeris::kepler(FLOAT M, FLOAT e)
{
    M = DEGREES_TO_RADIANS(M);
    
    FLOAT E = M;
    
    FLOAT previousE = E+1;
    
    for(int i=0; i<10; i++ )
    {
        E = E + (M+e*sin(E)-E)/(1-e*cos(E));
        
        // Optimize iterations
        if( previousE == E )
        {
            // No more iteration needed.
            break;
        }
        else
        {
            // Next iteration
            previousE = E;
        }
    }
    
    return RADIANS_TO_DEGREES(E);
}

#if !DISABLE_PLANETS
FLOAT Ephemeris::sumVSOP87Coefs(const VSOP87Coefficient *valuePlanetCoefficients, int coefCount, FLOAT T)
{
    // Parse each value in coef table
    FLOAT value = 0;
    for(int numCoef=0; numCoef<coefCount; numCoef++)
    {
        // Get coef
        VSOP87Coefficient coef;
        
#if ARDUINO
        // We limit SRAM usage by using flash memory (PROGMEM)
        memcpy_P(&coef, &valuePlanetCoefficients[numCoef], sizeof(VSOP87Coefficient));
#else
        coef = valuePlanetCoefficients[numCoef];
#endif
        
        FLOAT res = cos(coef.B + coef.C*T);
        
        // To avoid out of range issue with single precision
        // we've stored sqrt(A) and not A. As a result we need to square it back.
        res *= coef.A;
        res *= coef.A;
        
        value += res;
    }
    
    return value;
}
#endif

HorizontalCoordinates Ephemeris::equatorialToHorizontal(FLOAT H, FLOAT delta, FLOAT phi)
{
    HorizontalCoordinates coordinates;
    
    coordinates.azi = atan2(SIND(H), COSD(H)*SIND(phi)-TAND(delta)*COSD(phi));
    coordinates.azi = RADIANS_TO_DEGREES(coordinates.azi)+180; // +180 -> North is 0°
    coordinates.azi = LIMIT_DEGREES_TO_360(coordinates.azi);
    
    coordinates.alt = asin(SIND(phi)*SIND(delta) + COSD(phi)*COSD(delta)*COSD(H));
    coordinates.alt = RADIANS_TO_DEGREES(coordinates.alt);
    
    return coordinates;
}

EquatorialCoordinates Ephemeris::horizontalToEquatorial(FLOAT azimuth, FLOAT altitude, FLOAT latitude)
{
    EquatorialCoordinates coordinates;
    
    azimuth  = DEGREES_TO_RADIANS(azimuth-180); // -180 -> North is 0°
    altitude = DEGREES_TO_RADIANS(altitude);
    latitude = DEGREES_TO_RADIANS(latitude);
    
    
    coordinates.ra = atan2(sin(azimuth), cos(azimuth)*sin(latitude)+tan(altitude)*cos(latitude));
    coordinates.ra = RADIANS_TO_HOURS(coordinates.ra);
    coordinates.ra = LIMIT_HOURS_TO_24(coordinates.ra);
    
    coordinates.dec = asin(sin(latitude)*sin(altitude)-cos(latitude)*cos(altitude)*cos(azimuth));
    coordinates.dec = RADIANS_TO_DEGREES(coordinates.dec);
    
    return coordinates;
}

EquatorialCoordinates Ephemeris::EclipticToEquatorial(FLOAT lambda, FLOAT beta, FLOAT epsilon)
{
    lambda  = DEGREES_TO_RADIANS(lambda);
    beta    = DEGREES_TO_RADIANS(beta);
    epsilon = DEGREES_TO_RADIANS(epsilon);
    
    EquatorialCoordinates coordinates;
    coordinates.ra = atan2(sin(lambda)*cos(epsilon) - tan(beta)*sin(epsilon), cos(lambda));
    coordinates.ra = RADIANS_TO_HOURS(coordinates.ra);
    coordinates.ra = LIMIT_HOURS_TO_24(coordinates.ra);
    
    coordinates.dec = asin(sin(beta)*cos(epsilon) + cos(beta)*sin(epsilon)*sin(lambda));
    coordinates.dec = RADIANS_TO_DEGREES(coordinates.dec );
    
    return coordinates;
}

RectangularCoordinates Ephemeris::HeliocentricToRectangular(HeliocentricCoordinates hc, HeliocentricCoordinates hc0)
{
    RectangularCoordinates coordinates;
    
    coordinates.x = hc.radius * COSD(hc.lat) * COSD(hc.lon) - hc0.radius * COSD(hc0.lat) * COSD(hc0.lon);
    coordinates.y = hc.radius * COSD(hc.lat) * SIND(hc.lon) - hc0.radius * COSD(hc0.lat) * SIND(hc0.lon);
    coordinates.z = hc.radius * SIND(hc.lat)                - hc0.radius * SIND(hc0.lat);
    
    return coordinates;
}

FLOAT Ephemeris::meanGreenwichSiderealTimeAtDateAndTime(unsigned int day,   unsigned int month,   unsigned int year,
                                                        unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    FLOAT meanGreenwichSiderealTime;
    
    JulianDay jd0 = Calendar::julianDayForDateAndTime(day, month, year, 0, 0, 0);
    FLOAT T0        = T_WITH_JD(jd0.day,jd0.time);
    FLOAT T0Squared = T0*T0;
    FLOAT T0Cubed   = T0Squared*T0;
    
    // Sideral time at midnight
    FLOAT theta0 = 100.46061837 + (36000.770053608*T0) + (0.000387933*T0Squared) - (T0Cubed/38710000);
    theta0 = LIMIT_DEGREES_TO_360(theta0);
    theta0 = DEGREES_TO_HOURS(theta0);
    
    // Sideral time of day
    FLOAT thetaH = 1.00273790935*HOURS_MINUTES_SECONDS_TO_SECONDS(hours,minutes,seconds);
    thetaH = SECONDS_TO_DECIMAL_HOURS(thetaH);
    
    // Add time at midnight and time of day
    meanGreenwichSiderealTime = theta0 + thetaH;
    
    return LIMIT_HOURS_TO_24(meanGreenwichSiderealTime);
}

FLOAT Ephemeris::meanGreenwichSiderealTimeAtJD(JulianDay jd)
{
    FLOAT meanGreenwichSiderealTime;
    
    FLOAT T0        = T_WITH_JD(jd.day,jd.time);
    FLOAT T0Squared = T0*T0;
    FLOAT T0Cubed   = T0Squared*T0;
    
    // Sideral time at midnight
    FLOAT theta0 = 100.46061837 + (36000.770053608*T0) + (0.000387933*T0Squared) - (T0Cubed/38710000);
    theta0 = LIMIT_DEGREES_TO_360(theta0);
    theta0 = DEGREES_TO_HOURS(theta0);
    
    int day,month,year,hours,minutes,seconds;
    Calendar::dateAndTimeForJulianDay(jd, &day, &month, &year, &hours, &minutes, &seconds);
    
    // Sideral time of day
    FLOAT thetaH = 1.00273790935*HOURS_MINUTES_SECONDS_TO_SECONDS(hours,minutes,seconds);
    thetaH = SECONDS_TO_DECIMAL_HOURS(thetaH);
    
    // Add time at midnight and time of day
    meanGreenwichSiderealTime = theta0 + thetaH;
    
    return LIMIT_HOURS_TO_24(meanGreenwichSiderealTime);
}

#if !DISABLE_PLANETS
SolarSystemObject Ephemeris::solarSystemObjectAtDateAndTime(SolarSystemObjectIndex solarSystemObjectIndex,
                                                            unsigned int day,   unsigned int month,   unsigned int year,
                                                            unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    SolarSystemObject solarSystemObject;
    
    JulianDay jd = Calendar::julianDayForDateAndTime(day, month, year, hours, minutes, seconds);
    
    FLOAT T = T_WITH_JD(jd.day,jd.time);
    
    // Equatorial coordinates
    if( solarSystemObjectIndex == Sun )
    {
        solarSystemObject.equaCoordinates = equatorialCoordinatesForSunAtJD(jd, &solarSystemObject.distance);
    }
    else if( solarSystemObjectIndex == EarthsMoon )
    {
        solarSystemObject.equaCoordinates = equatorialCoordinatesForEarthsMoonAtJD(jd, &solarSystemObject.distance);
    }
    else
    {
        solarSystemObject.equaCoordinates = equatorialCoordinatesForPlanetAtJD(solarSystemObjectIndex, jd, &solarSystemObject.distance);
    }
    
    // Apparent diameter at a distance of 1 astronomical unit.
    FLOAT diameter = 0;
    switch (solarSystemObjectIndex)
    {
        case Mercury:
            diameter = 6.728;
            break;
            
        case Venus:
            diameter = 16.688;
            break;
            
        case Earth:
            diameter = NAN;
            break;
            
        case Mars:
            diameter = 9.364;
            break;
            
        case Jupiter:
            diameter = 197.146;
            break;
            
        case Saturn:
            diameter = 166.197;
            break;
            
        case Uranus:
            diameter = 70.476;
            break;
            
        case Neptune:
            diameter = 68.285;
            break;
            
        case Sun:
            diameter = 1919.26;
            break;
            
        case EarthsMoon:
            
            // 31.4'        -> 1884";
            // 385000.56 Km -> 385000.56*6.68459e-9 AU = 0.0025735709 AU
            // Theoretical diameter from a distance of 1 AU = 1884*0.0025735709 = 4.8486075756";
            diameter = 4.8486075756;
            
            break;
    }
    
    // Approximate apparent diameter in arc minutes according to distance
    solarSystemObject.diameter = diameter / solarSystemObject.distance/60;
    
    FLOAT meanSideralTime = meanGreenwichSiderealTimeAtDateAndTime(day, month, year, hours, minutes, seconds);
    
    FLOAT deltaNutation;
    FLOAT epsilon = obliquityAndNutationForT(T, NULL, &deltaNutation);
    
    // Apparent sideral time in floating hours
    FLOAT theta0 = meanSideralTime + (deltaNutation/15*COSD(epsilon))/3600;
    
    if( !isnan(longitudeOnEarth) && !isnan(latitudeOnEarth) )
    {
        // Geographic longitude in floating hours
        FLOAT L = DEGREES_TO_HOURS(longitudeOnEarth*longitudeOnEarthSign);
        
        // Geographic latitude in floating degrees
        FLOAT phi = latitudeOnEarth;
        
        // Local angle in floating degrees
        FLOAT H = (theta0-L-solarSystemObject.equaCoordinates.ra)*15;
        
        solarSystemObject.horiCoordinates = equatorialToHorizontal(H,solarSystemObject.equaCoordinates.dec,phi);
        
        // Mean sideral time at midnight
        FLOAT T0 = Ephemeris::meanGreenwichSiderealTimeAtDateAndTime(day, month, year, 0, 0, 0);
        
        // Compute rise and set
        EquatorialCoordinates tmpCoord,tmpCoord0,tmpCoord24;
        FLOAT linearSpeedRA,linearSpeedDec;
        switch (solarSystemObjectIndex)
        {
            case Sun:
                
                //
                // Assume Sun speed to be linear for 24 hour range
                //
                
                tmpCoord0  = equatorialCoordinatesForSunAtJD(Calendar::julianDayForDateAndTime(day, month, year,  0, 0, 0), NULL);
                tmpCoord24 = equatorialCoordinatesForSunAtJD(Calendar::julianDayForDateAndTime(day, month, year, 24, 0, 0), NULL);
                
                linearSpeedRA  = (tmpCoord24.ra  - tmpCoord0.ra);
                linearSpeedDec = (tmpCoord24.dec - tmpCoord0.dec);
                
                // First approximation at midday
                tmpCoord.ra  = tmpCoord0.ra  + linearSpeedRA*0.5;
                tmpCoord.dec = tmpCoord0.dec + linearSpeedDec*0.5;
                solarSystemObject.riseAndSetState = riseAndSetForEquatorialCoordinatesAndT0(tmpCoord,
                                                                                            T0,
                                                                                            &solarSystemObject.rise, &solarSystemObject.set,
                                                                                            0,
                                                                                            solarSystemObject.diameter/60.0);
                if( !isnan(solarSystemObject.rise) )
                {
                    // Now interpolate coordinates at rise time estimation
                    tmpCoord.ra  = tmpCoord0.ra  + linearSpeedRA*solarSystemObject.rise/24.0;
                    tmpCoord.dec = tmpCoord0.dec + linearSpeedDec*solarSystemObject.rise/24.0;
                    
                    // Compute new coordinates to improve precision
                    riseAndSetForEquatorialCoordinatesAndT0(tmpCoord,
                                                            T0,
                                                            &solarSystemObject.rise, NULL,
                                                            0,
                                                            solarSystemObject.diameter/60.0);
                }
                
                if( !isnan(solarSystemObject.set) )
                {
                    // Now interpolate coordinates at set time estimation
                    tmpCoord.ra  = tmpCoord0.ra  + linearSpeedRA*solarSystemObject.set/24.0;
                    tmpCoord.dec = tmpCoord0.dec + linearSpeedDec*solarSystemObject.set/24.0;
                    
                    // Compute new coordinates to improve precision
                    riseAndSetForEquatorialCoordinatesAndT0(tmpCoord,
                                                            T0,
                                                            NULL, &solarSystemObject.set,
                                                            0,
                                                            solarSystemObject.diameter/60.0);
                }
                
                break;
                
            case EarthsMoon:

                solarSystemObject.riseAndSetState = riseAndSetForEquatorialCoordinatesAndT0(solarSystemObject.equaCoordinates,
                                                                                            T0,
                                                                                            &solarSystemObject.rise, &solarSystemObject.set,
                                                                                            0,
                                                                                            0);
                
                //
                // Assume Moon speed to be linear for 24 hour range
                //
                
                tmpCoord0  = equatorialCoordinatesForEarthsMoonAtJD(Calendar::julianDayForDateAndTime(day, month, year, hours-12, minutes, seconds), NULL);
                tmpCoord24 = equatorialCoordinatesForEarthsMoonAtJD(Calendar::julianDayForDateAndTime(day, month, year, hours+12, minutes, seconds), NULL);
                
                linearSpeedRA  = (tmpCoord24.ra  - tmpCoord0.ra)/24.f;
                linearSpeedDec = (tmpCoord24.dec - tmpCoord0.dec)/24.f;
                
                // First approximation at midday
                solarSystemObject.riseAndSetState = riseAndSetForEquatorialCoordinatesAndT0(solarSystemObject.equaCoordinates,
                                                                                            T0,
                                                                                            &solarSystemObject.rise, &solarSystemObject.set,
                                                                                            0,
                                                                                            0);
                
                if( !isnan(solarSystemObject.rise) )
                {
                    // Now interpolate coordinates at rise time estimation
                    tmpCoord.ra  = solarSystemObject.equaCoordinates.ra  + linearSpeedRA*solarSystemObject.rise;
                    tmpCoord.dec = solarSystemObject.equaCoordinates.dec + linearSpeedDec*solarSystemObject.rise;
                    
                    // Compute new coordinates to improve precision
                    riseAndSetForEquatorialCoordinatesAndT0(tmpCoord,
                                                            T0,
                                                            &solarSystemObject.rise, NULL,
                                                            57/60.0,
                                                            solarSystemObject.diameter/60.0);
                }
                
                if( !isnan(solarSystemObject.set) )
                {
                    // Now interpolate coordinates at set time estimation
                    tmpCoord.ra  = solarSystemObject.equaCoordinates.ra  + linearSpeedRA*solarSystemObject.set;
                    tmpCoord.dec = solarSystemObject.equaCoordinates.dec + linearSpeedDec*solarSystemObject.set;
                    
                    // Compute new coordinates to improve precision
                    riseAndSetForEquatorialCoordinatesAndT0(tmpCoord,
                                                            T0,
                                                            NULL, &solarSystemObject.set,
                                                            57/60.0,
                                                            solarSystemObject.diameter/60.0);
                }
                
                break;
                
            default:
                solarSystemObject.riseAndSetState = riseAndSetForEquatorialCoordinatesAndT0(solarSystemObject.equaCoordinates,
                                                                                            T0,
                                                                                            &solarSystemObject.rise, &solarSystemObject.set,
                                                                                            0,
                                                                                            0);
                break;
        }
    }
    else
    {
        solarSystemObject.horiCoordinates.alt = NAN;
        solarSystemObject.horiCoordinates.azi = NAN;
    }
    
    return solarSystemObject;
}
#else
SolarSystemObject Ephemeris::solarSystemObjectAtDateAndTime(SolarSystemObjectIndex solarSystemObjectIndex,
                                                            unsigned int day,   unsigned int month,   unsigned int year,
                                                            unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    // If DISABLE_PLANETS we simply return an empty SolarSystemObject
    
    SolarSystemObject solarSystemObject;
    memset(&solarSystemObject,'\0',sizeof(SolarSystemObject));
    return solarSystemObject;
}
#endif

void Ephemeris::setLocationOnEarth(FLOAT floatingLatitude, FLOAT floatingLongitude)
{
    latitudeOnEarth  = floatingLatitude;
    longitudeOnEarth = floatingLongitude;
}

void Ephemeris::setLocationOnEarth(FLOAT latDegrees, FLOAT latMinutes, FLOAT latSeconds,
                                   FLOAT lonDegrees, FLOAT lonMinutes, FLOAT lonSeconds)
{
    latitudeOnEarth  = DEGREES_MINUTES_SECONDS_TO_DECIMAL_DEGREES(latDegrees,latMinutes,latSeconds);
    longitudeOnEarth = DEGREES_MINUTES_SECONDS_TO_DECIMAL_DEGREES(lonDegrees,lonMinutes,lonSeconds);
}

#if !DISABLE_PLANETS
EquatorialCoordinates  Ephemeris::equatorialCoordinatesForPlanetAtJD(SolarSystemObjectIndex solarSystemObjectIndex, JulianDay jd, FLOAT *distance)
{
    EquatorialCoordinates coordinates;
    coordinates.ra  = 0;
    coordinates.dec = 0;
    
    FLOAT T      = T_WITH_JD(jd.day,jd.time);
    FLOAT lastT  = 0;
    FLOAT TLight = 0;
    HeliocentricCoordinates hcPlanet;
    HeliocentricCoordinates hcEarth;
    RectangularCoordinates  rectPlanet;
    
    JulianDay jd2;
    jd2.day  = jd.day;
    jd2.time = jd.time;
    
    FLOAT x2=0,y2=0,z2=0;
    
    FLOAT dist = 0;
    
    // Iterate for good precision according to light speed delay
    while (T != lastT)
    {
        lastT = T;
        
        FLOAT jd2 = (jd.day+jd.time) - TLight;
        T = T_WITH_JD(jd2,0);
        
        hcPlanet   = Ephemeris::heliocentricCoordinatesForPlanetAndT(solarSystemObjectIndex, T);
        if( isnan(hcPlanet.radius)  )
        {
            break;
        }
        
        hcEarth    = Ephemeris::heliocentricCoordinatesForPlanetAndT(Earth, T);
        
        rectPlanet = HeliocentricToRectangular(hcPlanet,hcEarth);
        
        // Precomputed square
        x2 = rectPlanet.x*rectPlanet.x;
        y2 = rectPlanet.y*rectPlanet.y;
        z2 = rectPlanet.z*rectPlanet.z;
        
        // Real distance from Earth
        dist = sqrtf(x2+y2+z2);

        // Light time
        TLight = dist * 0.0057755183;
    }
    
    if( distance ) *distance = dist;

    // Geocentic longitude
    FLOAT lambda = atan2(rectPlanet.y,rectPlanet.x);
    lambda = RADIANS_TO_DEGREES(lambda);
    lambda = LIMIT_DEGREES_TO_360(lambda);
    
    // Geocentric latitude
    FLOAT beta   = rectPlanet.z/sqrtf(x2+y2);
    beta = RADIANS_TO_DEGREES(beta);
    
    
    // Remove abberation
    {
        PlanetayOrbit earthOrbit = Ephemeris::planetayOrbitForPlanetAndT(Earth, T);
        
        FLOAT TSquared = T*T;
        
        FLOAT L0 = 280.46646 + T*36000.76983 + TSquared*0.0003032;
        L0 = LIMIT_DEGREES_TO_360(L0);
        
        FLOAT M = 357.52911 + T*35999.05029  - TSquared*0.0001537;
        M = LIMIT_DEGREES_TO_360(M);
        
        FLOAT C =
        +(1.914602 - T*0.004817 - TSquared*0.000014) * SIND(M)
        +(0.019993 - T*0.000101                    ) * SIND(2*M)
        + 0.000289                                   * SIND(3*M);
        
        // Sun longitude
        FLOAT O = L0 + C;
        
        // Abberation
        FLOAT k = 20.49552;
        FLOAT xAberration = (-k*COSD(O - lambda) + earthOrbit.e*k*COSD(earthOrbit.pi - lambda)) / COSD(beta)/3600;
        FLOAT yAberration = -k*SIND(beta)*(SIND(O - lambda) - earthOrbit.e*SIND(earthOrbit.pi - lambda))/3600;
        lambda -= xAberration;
        beta   -= yAberration;
    }
    
    // Obliquity and Nutation
    FLOAT deltaNutation;
    FLOAT epsilon = Ephemeris::obliquityAndNutationForT(T, NULL, &deltaNutation);
    
    // Intergrate nutation
    lambda += deltaNutation/3600;
    
    EquatorialCoordinates eqCoord = EclipticToEquatorial(lambda, beta, epsilon);
    
    //
    // Geocentric to topocentric conversion
    //
    
    FLOAT paralax = 8.794/(dist)/3600.f;
    
    FLOAT meanSideralTime = meanGreenwichSiderealTimeAtJD(jd);
    
    // Geographic longitude in floating hours
    FLOAT L = DEGREES_TO_HOURS(longitudeOnEarth*longitudeOnEarthSign);
    
    // Apparent sideral time in floating hours
    FLOAT theta0 = meanSideralTime + (deltaNutation/15*COSD(epsilon))/3600;
    
    // Local angle in floating degrees
    FLOAT H = (theta0-L-eqCoord.ra)*15;
    
    // Geocentric rectangular coordinates of the observer (Chap 6 p33)
    FLOAT u = ATAND(0.99664719*TAND(latitudeOnEarth));
    FLOAT rhoSinDeltaPrime    = 0.99664719*SIND(u)+(altitudeOnEarth/6378140)*SIND(latitudeOnEarth);
    FLOAT rhoCosDeltaPrime = COSD(u)+(altitudeOnEarth/6378140)*COSD(latitudeOnEarth);
    
    // Paralax correction (Chap 27 p 105)
    FLOAT detlaRa = ATAND(-rhoCosDeltaPrime*SIND(paralax)*SIND(H)/(COSD(eqCoord.dec)-rhoCosDeltaPrime*SIND(paralax)*COSD(H)));
    FLOAT detlaDec = ATAND((SIND(eqCoord.dec)-rhoSinDeltaPrime*SIND(paralax))*COSD(detlaRa)/(COSD(eqCoord.dec)-rhoCosDeltaPrime*SIND(paralax)*COSD(H)));
    
    eqCoord.ra += detlaRa/15.f;
    eqCoord.dec = detlaDec;
    
    return eqCoord;
}
#endif

#if !DISABLE_PLANETS
HeliocentricCoordinates  Ephemeris::heliocentricCoordinatesForPlanetAndT(SolarSystemObjectIndex solarSystemObjectIndex, FLOAT T)
{
    HeliocentricCoordinates coordinates;
    
    T = T/10;
    FLOAT TSquared = T*T;
    FLOAT TCubed   = TSquared*T;
    FLOAT T4       = TCubed*T;
    FLOAT T5       = T4*T;
    
    int SizeOfVSOP87Coefficient = sizeof(VSOP87Coefficient);
    
    FLOAT l0=0,l1=0,l2=0,l3=0,l4=0,l5=0;
    FLOAT b0=0,b1=0,b2=0,b3=0,b4=0,b5=0;
    FLOAT r0=0,r1=0,r2=0,r3=0,r4=0,r5=0;
    
    switch (solarSystemObjectIndex)
    {
        case Sun:
            break;
            
        case Mercury:
            
            l0 = sumVSOP87Coefs(L0MercuryCoefficients,sizeof(L0MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1MercuryCoefficients,sizeof(L1MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2MercuryCoefficients,sizeof(L2MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3MercuryCoefficients,sizeof(L3MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4MercuryCoefficients,sizeof(L4MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            l5 = sumVSOP87Coefs(L5MercuryCoefficients,sizeof(L5MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0MercuryCoefficients,sizeof(B0MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1MercuryCoefficients,sizeof(B1MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            b2 = sumVSOP87Coefs(B2MercuryCoefficients,sizeof(B2MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            b3 = sumVSOP87Coefs(B3MercuryCoefficients,sizeof(B3MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            b4 = sumVSOP87Coefs(B4MercuryCoefficients,sizeof(B4MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0MercuryCoefficients,sizeof(R0MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1MercuryCoefficients,sizeof(R1MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2MercuryCoefficients,sizeof(R2MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3MercuryCoefficients,sizeof(R3MercuryCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
            
        case Venus:
            
            l0 = sumVSOP87Coefs(L0VenusCoefficients,sizeof(L0VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1VenusCoefficients,sizeof(L1VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2VenusCoefficients,sizeof(L2VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3VenusCoefficients,sizeof(L3VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4VenusCoefficients,sizeof(L4VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            l5 = sumVSOP87Coefs(L5VenusCoefficients,sizeof(L5VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0VenusCoefficients,sizeof(B0VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1VenusCoefficients,sizeof(B1VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            b2 = sumVSOP87Coefs(B2VenusCoefficients,sizeof(B2VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            b3 = sumVSOP87Coefs(B3VenusCoefficients,sizeof(B3VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            b4 = sumVSOP87Coefs(B4VenusCoefficients,sizeof(B4VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0VenusCoefficients,sizeof(R0VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1VenusCoefficients,sizeof(R1VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2VenusCoefficients,sizeof(R2VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3VenusCoefficients,sizeof(R3VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            r4 = sumVSOP87Coefs(R4VenusCoefficients,sizeof(R4VenusCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
            
        case Earth:
            
            l0 = sumVSOP87Coefs(L0EarthCoefficients,sizeof(L0EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1EarthCoefficients,sizeof(L1EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2EarthCoefficients,sizeof(L2EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3EarthCoefficients,sizeof(L3EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4EarthCoefficients,sizeof(L4EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            l5 = sumVSOP87Coefs(L5EarthCoefficients,sizeof(L5EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0EarthCoefficients,sizeof(B0EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1EarthCoefficients,sizeof(B1EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0EarthCoefficients,sizeof(R0EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1EarthCoefficients,sizeof(R1EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2EarthCoefficients,sizeof(R2EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3EarthCoefficients,sizeof(R3EarthCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
            
        case Mars:
            
            l0 = sumVSOP87Coefs(L0MarsCoefficients,sizeof(L0MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1MarsCoefficients,sizeof(L1MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2MarsCoefficients,sizeof(L2MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3MarsCoefficients,sizeof(L3MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4MarsCoefficients,sizeof(L4MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            l5 = sumVSOP87Coefs(L5MarsCoefficients,sizeof(L5MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0MarsCoefficients,sizeof(B0MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1MarsCoefficients,sizeof(B1MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            b2 = sumVSOP87Coefs(B2MarsCoefficients,sizeof(B2MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            b3 = sumVSOP87Coefs(B3MarsCoefficients,sizeof(B3MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            b4 = sumVSOP87Coefs(B4MarsCoefficients,sizeof(B4MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0MarsCoefficients,sizeof(R0MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1MarsCoefficients,sizeof(R1MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2MarsCoefficients,sizeof(R2MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3MarsCoefficients,sizeof(R3MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            r4 = sumVSOP87Coefs(R4MarsCoefficients,sizeof(R4MarsCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
            
        case Jupiter:
            
            l0 = sumVSOP87Coefs(L0JupiterCoefficients,sizeof(L0JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1JupiterCoefficients,sizeof(L1JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2JupiterCoefficients,sizeof(L2JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3JupiterCoefficients,sizeof(L3JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4JupiterCoefficients,sizeof(L4JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            l5 = sumVSOP87Coefs(L5JupiterCoefficients,sizeof(L5JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0JupiterCoefficients,sizeof(B0JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1JupiterCoefficients,sizeof(B1JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            b2 = sumVSOP87Coefs(B2JupiterCoefficients,sizeof(B2JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            b3 = sumVSOP87Coefs(B3JupiterCoefficients,sizeof(B3JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            b4 = sumVSOP87Coefs(B4JupiterCoefficients,sizeof(B4JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            b5 = sumVSOP87Coefs(B5JupiterCoefficients,sizeof(B5JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0JupiterCoefficients,sizeof(R0JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1JupiterCoefficients,sizeof(R1JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2JupiterCoefficients,sizeof(R2JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3JupiterCoefficients,sizeof(R3JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            r4 = sumVSOP87Coefs(R4JupiterCoefficients,sizeof(R4JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            r5 = sumVSOP87Coefs(R5JupiterCoefficients,sizeof(R5JupiterCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
            
        case Saturn:
            
            l0 = sumVSOP87Coefs(L0SaturnCoefficients,sizeof(L0SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1SaturnCoefficients,sizeof(L1SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2SaturnCoefficients,sizeof(L2SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3SaturnCoefficients,sizeof(L3SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4SaturnCoefficients,sizeof(L4SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            l5 = sumVSOP87Coefs(L5SaturnCoefficients,sizeof(L5SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0SaturnCoefficients,sizeof(B0SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1SaturnCoefficients,sizeof(B1SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            b2 = sumVSOP87Coefs(B2SaturnCoefficients,sizeof(B2SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            b3 = sumVSOP87Coefs(B3SaturnCoefficients,sizeof(B3SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            b4 = sumVSOP87Coefs(B4SaturnCoefficients,sizeof(B4SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            b5 = sumVSOP87Coefs(B5SaturnCoefficients,sizeof(B5SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0SaturnCoefficients,sizeof(R0SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1SaturnCoefficients,sizeof(R1SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2SaturnCoefficients,sizeof(R2SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3SaturnCoefficients,sizeof(R3SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            r4 = sumVSOP87Coefs(R4SaturnCoefficients,sizeof(R4SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            r5 = sumVSOP87Coefs(R5SaturnCoefficients,sizeof(R5SaturnCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
            
        case Uranus:
            
            l0 = sumVSOP87Coefs(L0UranusCoefficients,sizeof(L0UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1UranusCoefficients,sizeof(L1UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2UranusCoefficients,sizeof(L2UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3UranusCoefficients,sizeof(L3UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4UranusCoefficients,sizeof(L4UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0UranusCoefficients,sizeof(B0UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1UranusCoefficients,sizeof(B1UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            b2 = sumVSOP87Coefs(B2UranusCoefficients,sizeof(B2UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            b3 = sumVSOP87Coefs(B3UranusCoefficients,sizeof(B3UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            b4 = sumVSOP87Coefs(B4UranusCoefficients,sizeof(B4UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0UranusCoefficients,sizeof(R0UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1UranusCoefficients,sizeof(R1UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2UranusCoefficients,sizeof(R2UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3UranusCoefficients,sizeof(R3UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            r4 = sumVSOP87Coefs(R4UranusCoefficients,sizeof(R4UranusCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
            
        case Neptune:
            
            l0 = sumVSOP87Coefs(L0NeptuneCoefficients,sizeof(L0NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            l1 = sumVSOP87Coefs(L1NeptuneCoefficients,sizeof(L1NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            l2 = sumVSOP87Coefs(L2NeptuneCoefficients,sizeof(L2NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            l3 = sumVSOP87Coefs(L3NeptuneCoefficients,sizeof(L3NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            l4 = sumVSOP87Coefs(L4NeptuneCoefficients,sizeof(L4NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            
            b0 = sumVSOP87Coefs(B0NeptuneCoefficients,sizeof(B0NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            b1 = sumVSOP87Coefs(B1NeptuneCoefficients,sizeof(B1NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            b2 = sumVSOP87Coefs(B2NeptuneCoefficients,sizeof(B2NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            b3 = sumVSOP87Coefs(B3NeptuneCoefficients,sizeof(B3NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            b4 = sumVSOP87Coefs(B4NeptuneCoefficients,sizeof(B4NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            
            r0 = sumVSOP87Coefs(R0NeptuneCoefficients,sizeof(R0NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            r1 = sumVSOP87Coefs(R1NeptuneCoefficients,sizeof(R1NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            r2 = sumVSOP87Coefs(R2NeptuneCoefficients,sizeof(R2NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            r3 = sumVSOP87Coefs(R3NeptuneCoefficients,sizeof(R3NeptuneCoefficients)/SizeOfVSOP87Coefficient,T);
            
            break;
            
        default:
            
            // Do not work for Moon...
            l0 = NAN;
            l1 = NAN;
            l2 = NAN;
            l3 = NAN;
            l4 = NAN;
            l5 = NAN;
            
            b0 = NAN;
            b1 = NAN;
            b2 = NAN;
            b3 = NAN;
            b4 = NAN;
            b5 = NAN;
            
            r0 = NAN;
            r1 = NAN;
            r2 = NAN;
            r3 = NAN;
            r4 = NAN;
            r5 = NAN;
            
            break;
    }
    
    // L
    coordinates.lon  = (l0 + l1*T + l2*TSquared + l3*TCubed + l4*T4 + l5*T5)/100000000.0;
    coordinates.lon  = RADIANS_TO_DEGREES(coordinates.lon);
    coordinates.lon  = LIMIT_DEGREES_TO_360(coordinates.lon);
    
    // B
    coordinates.lat  = (b0 + b1*T + b2*TSquared + b3*TCubed + b4*T4 + b5*T5)/100000000.0;
    coordinates.lat  = RADIANS_TO_DEGREES(coordinates.lat);
    
    // R
    coordinates.radius = (r0 + r1*T + r2*TSquared + r3*TCubed + r4*T4 + r5*T5)/100000000.0;
    
    return coordinates;
}
#endif

RiseAndSetState  Ephemeris::riseAndSetForEquatorialCoordinatesAndT0(EquatorialCoordinates coord, FLOAT T0, FLOAT *rise, FLOAT *set,
                                                                    FLOAT paralax, FLOAT apparentDiameter)
{
    if( rise ) *rise = NAN;
    if( set )  *set  = NAN;
    
    if( isnan(longitudeOnEarth) && isnan(latitudeOnEarth) )
    {
        return LocationOnEarthUnitialized;
    }
    
    FLOAT lon = DEGREES_TO_HOURS(longitudeOnEarth*longitudeOnEarthSign);
    FLOAT lat = latitudeOnEarth;
    
    // Altitude angle
    FLOAT n1 = ACOSD(6378140/(FLOAT)(6378140 + altitudeOnEarth));
    
    // h0 = P - R - 1/2 d - η1
    // P: Parallax in degrees (57' for the moon and 0 for others).
    // R: Refraction of Radau in degrees
    // d: Apparent diameter in degrees (32' for moon en sun).
    // n1 Altitude parameter
    FLOAT h0 = paralax - 34/60.0 -apparentDiameter/2 - n1;
    
    // cos H = (sin(h0) - sin(φ) sin(δ))/(cos(φ) cos (δ))
    FLOAT cosH = (SIND(h0) - SIND(lat)*SIND(coord.dec))/(COSD(lat)*COSD(coord.dec));
    
    if( cosH < -1 )
    {
        return ObjectAlwaysInSky;
    }
    else if( cosH > 1 )
    {
        return ObjectNeverInSky;
    }
    
    FLOAT H = ACOSD(cosH);
    H = DEGREES_TO_HOURS(H);
    
    if( rise )
    {
        *rise = coord.ra - H + lon - T0;
        
        // Sideral time to mean time
        *rise *= 0.9972695663;
    }
    
    if( set )
    {
        *set  = coord.ra + H + lon - T0;
        
        // Sideral time to mean time
        *set  *= 0.9972695663;
    }
    
    return RiseAndSetOk;
}

RiseAndSetState Ephemeris::riseAndSetForEquatorialCoordinatesAtDateAndTime(EquatorialCoordinates coord,
                                                                           FLOAT *rise, FLOAT *set,
                                                                           unsigned int day,  unsigned int month,  unsigned int year,
                                                                           unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    // Mean sideral time at midnight
    FLOAT T0 = Ephemeris::meanGreenwichSiderealTimeAtDateAndTime(day, month, year, 0, 0, 0);
    
    return riseAndSetForEquatorialCoordinatesAndT0(coord, T0, rise, set, 0, 0);
}

void Ephemeris::setAltitude(int altitude)
{
    altitudeOnEarth = altitude;
}

void Ephemeris::flipLongitude(bool flip)
{
    if( flip == true )
    {
        longitudeOnEarthSign = 1;
    }
    else
    {
        longitudeOnEarthSign = -1;
    }
}

FLOAT Ephemeris::floatingHoursWithUTCOffset(float floatingHours, int UTCOffset)
{
    floatingHours += UTCOffset;
    
    return LIMIT_HOURS_TO_24(floatingHours);
}

EquatorialCoordinates Ephemeris::equatorialEquinoxToEquatorialJNowAtDateAndTime(EquatorialCoordinates eqEquinoxCoordinates,
                                                                                int equinox,
                                                                                unsigned int day,  unsigned int month,  unsigned int year,
                                                                                unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    EquatorialCoordinates noDrift;
    noDrift.ra  = 0;
    noDrift.dec = 0;
    
    return equatorialEquinoxToEquatorialJNowAtDateAndTime(eqEquinoxCoordinates,
                                                          equinox,
                                                          noDrift,
                                                          day,month,year,
                                                          hours,minutes,seconds);
}

EquatorialCoordinates Ephemeris::equatorialEquinoxToEquatorialJNowAtDateAndTime(EquatorialCoordinates eqEquinoxCoordinates,
                                                                                int equinox,
                                                                                EquatorialCoordinates eqDriftPerYear,
                                                                                unsigned int day,   unsigned int month,   unsigned int year,
                                                                                unsigned int hours, unsigned int minutes, unsigned int seconds)
{
    JulianDay jd = Calendar::julianDayForDateAndTime(day, month, year, hours, minutes, seconds);
    
    FLOAT T = T_WITH_JD(jd.day,jd.time);
    
    return equatorialEquinoxToEquatorialJNowAtDateForT(eqEquinoxCoordinates, equinox, eqDriftPerYear, T, year);
}

EquatorialCoordinates Ephemeris::equatorialEquinoxToEquatorialJNowAtDateForT(EquatorialCoordinates eqEquinoxCoordinates,
                                                                            int equinox,
                                                                            EquatorialCoordinates eqDriftPerYear,
                                                                            FLOAT T,
                                                                            unsigned int year)
{
    
    EquatorialCoordinates eqJNowCoordinates = eqEquinoxCoordinates;
    
    
    
    //
    // Precession of the equinoxes
    //
    
    // Convert RA to degrees
    FLOAT RADeg = HOURS_TO_DEGREES(eqEquinoxCoordinates.ra);
    
#if 1
    FLOAT Teq  = equinox-2000;
    FLOAT m    = 3.07496 + 0.00186*Teq; // S
    FLOAT nRA  = 1.33621 - 0.00057*Teq; // S
    FLOAT nDec = 20.0431 - 0.0085*Teq;  // "
    
    FLOAT deltaRA  = m + nRA * SIND(RADeg)*TAND(eqEquinoxCoordinates.dec);
    FLOAT deltaDec = nDec * COSD(RADeg);
#else

    FLOAT deltaRA  = 0;
    FLOAT deltaDec = 0;
    
    // TODO: more precise method.
    
#endif
    
    // Add self movement per year
    deltaRA  += eqDriftPerYear.ra;
    deltaDec += eqDriftPerYear.dec;
    
    eqJNowCoordinates.ra  += SECONDS_TO_DECIMAL_HOURS(deltaRA*(year-equinox));
    eqJNowCoordinates.dec += SECONDS_TO_DECIMAL_DEGREES(deltaDec*(year-equinox));
    
    
    //
    // Nutation
    //
    
    // For tests
    //eqJNowCoordinates.ra  = hoursMinutesSecondsToFloatingHours(10,9,7.845);
    //eqJNowCoordinates.dec = degreesMinutesSecondsToFloatingDegrees(99, 53, 48.96);
    
    RADeg  = HOURS_TO_DEGREES(eqJNowCoordinates.ra);
    FLOAT dec = eqJNowCoordinates.dec;
    
    FLOAT deltaPhi; // Delta nutation
    FLOAT deltaEps; // Delta obliquity
    FLOAT eps = obliquityAndNutationForT(T, &deltaEps, &deltaPhi);
    
    FLOAT deltaNutationRA = (COSD(eps) + SIND(eps) * SIND(RADeg) * TAND(dec))*deltaPhi - (COSD(RADeg)*TAND(dec)*deltaEps);
    deltaNutationRA = SECONDS_TO_DECIMAL_DEGREES(deltaNutationRA);
    deltaNutationRA = DEGREES_TO_HOURS(deltaNutationRA);
    
    FLOAT deltaNutationDec = (SIND(eps)*COSD(RADeg))*deltaPhi + SIND(RADeg)*deltaEps;
    deltaNutationDec = SECONDS_TO_DECIMAL_DEGREES(deltaNutationDec);
    
    eqJNowCoordinates.ra  += deltaNutationRA;
    eqJNowCoordinates.dec += deltaNutationDec;
    
    
    //
    // Aberration
    //
    
    FLOAT TSquared = T*T;
    
    FLOAT L0 = 280.46646 + T*36000.76983 + TSquared*0.0003032;
    L0 = LIMIT_DEGREES_TO_360(L0);
    
    FLOAT M = 357.5291092 + T*35999.0502909  - TSquared*0.0001536;
    M = LIMIT_DEGREES_TO_360(M);
    
    FLOAT e = 0.016708634 - T*0.000042037 - TSquared*0.0000001267;
    
    FLOAT pi = 102.93735 + T*1.71946 + TSquared*0.00046;
    
    FLOAT C =
    +(1.914602 - T*0.004817 - TSquared*0.000014) * SIND(M)
    +(0.019993 - T*0.000101                    ) * SIND(2*M)
    + 0.000289                                   * SIND(3*M);
    
    FLOAT O = L0 + C;
    
    FLOAT K = 20.49552;
    
    RADeg  = HOURS_TO_DEGREES(eqJNowCoordinates.ra);
    dec    = eqJNowCoordinates.dec;
    
    FLOAT deltaAberrationRA =
    -K  *(COSD(RADeg)*COSD(O)*COSD(eps)+SIND(RADeg)*SIND(O))/COSD(dec)
    +K*e*(COSD(RADeg)*COSD(pi)*COSD(eps)+SIND(RADeg)*SIND(pi))/COSD(dec);
    deltaAberrationRA = SECONDS_TO_DECIMAL_DEGREES(deltaAberrationRA);
    deltaAberrationRA = DEGREES_TO_HOURS(deltaAberrationRA);
    
    FLOAT deltaAberrationDec =
    -K  *(COSD(O)*COSD(eps)*(TAND(eps)*COSD(dec)-SIND(RADeg)*SIND(dec)) + COSD(RADeg)*SIND(dec)*SIND(O))
    +K*e*(COSD(pi)*COSD(eps)*(TAND(eps)*COSD(dec)-SIND(RADeg)*SIND(dec))+ COSD(RADeg)*SIND(dec)*SIND(pi));
    deltaAberrationDec = SECONDS_TO_DECIMAL_DEGREES(deltaAberrationDec);
    
    eqJNowCoordinates.ra  += deltaAberrationRA;
    eqJNowCoordinates.dec += deltaAberrationDec;
    
    
    //
    // Avoid overflow
    //
    
    if( eqJNowCoordinates.dec>90 )
    {
        eqJNowCoordinates.dec = 180-eqJNowCoordinates.dec;
        eqJNowCoordinates.ra += 12;
    }
    else if( eqJNowCoordinates.dec<-90 )
    {
        eqJNowCoordinates.dec = -180-eqJNowCoordinates.dec;
        eqJNowCoordinates.ra += 12;
    }
    
    eqJNowCoordinates.ra  = LIMIT_HOURS_TO_24(eqJNowCoordinates.ra);
    
    return eqJNowCoordinates;
}
