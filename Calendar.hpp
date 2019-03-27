/*
 * Calendar.hpp
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

#ifndef Calendar_h
#define Calendar_h

#ifndef FLOAT
// Change FLOAT according to wanted precision (float or double).
#define FLOAT float
#endif


/*! This structure describes a Julian day with integer and floating value parts.
 *  As a result we can get 64bits precision with single precision processor. */
struct JulianDay
{
    /*! Integer value of Jullian day. */
    long day;
    
    /*! Floating value of Jullian day (0.0=0h 1.0=24h) */
    FLOAT time;
};


/*!
 * This class manipulates calendar arithmetic.
 */
class Calendar
{
public:
    
    static JulianDay julianDayForDate(FLOAT day, int month, int year);
    
    static JulianDay julianDayForDateAndTime(int day,   int month,   int year,
                                             int hours, int minutes, int seconds);
    
    static void dateForJulianDay(JulianDay julianDay, FLOAT *day, int *month, int *year);
    
    static void dateAndTimeForJulianDay(JulianDay julianDay,
                                        int *day,   int *month,   int *year,
                                        int *hours, int *minutes, int *seconds);
    
    static unsigned int weekDayForDate( int day, int month, int year);
    
    static unsigned int weekDayForJulianDay(JulianDay julianDay);
};

#endif
