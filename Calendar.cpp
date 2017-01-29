/*
 * Calendar.cpp
 */

#if ARDUINO
#include <Arduino.h>
#endif

#include "Calendar.hpp"

#define INT(value) (value >= 0 ? (long)(value) : (long)(value-1))

JulianDay Calendar::julianDayForDate(float day, unsigned int month, unsigned int year )
{
    if( month <= 2 )
    {
        year = year - 1;
        month = month + 12;
    }
    
    long A = INT(year/100.0);
    long B = 2-A+INT(A/4);
    
    JulianDay julianDay;
    
    julianDay.day  = INT(365.25*(year+4716)) + INT(30.6001*(month+1)) + (int)(day) + B - 1524.5;
    
    julianDay.time = 0.5 + (day-(int)day);
    
    // Adjust day and time if needed
    if( julianDay.time >= 1 )
    {
        julianDay.time -=1;
        
        julianDay.day += 1;
    }
    
    return julianDay;
}

JulianDay Calendar::julianDayForDateAndTime(unsigned int day, unsigned int month, unsigned int year,
                                            unsigned int hour, unsigned int minute, unsigned int second)
{
    return julianDayForDate(day+ hour/24.0 + minute/1440.0 + second/86400.0,month,year);
}


void Calendar::dateForJulianDay(JulianDay julianDay, float *day, unsigned int *month, unsigned int *year )
{
    long Z  = julianDay.day;
    
    float F = julianDay.time + 0.5;
    
    // Adjust day and time if needed
    if( F >= 1 )
    {
        F -= 1;
        Z += 1;
    }
    
    long A;
    if( Z<2299161 )
    {
        A = Z;
    }
    else
    {
        long alpha = INT((Z-1867216.25)/36524.25);
        A     = Z+1+alpha-INT(alpha/4);
    }
    
    long B = A + 1524;
    long C = INT((B-122.1)/365.25);
    long D = INT(365.25*C);
    long E = INT((B-D)/30.6001);
    
    *day   = B-D-INT(30.6001*E)+F;
    if( E<14)
    {
        *month = (unsigned)(E-1);
    }
    else
    {
        *month = (unsigned)(E-13);
    }
    
    if( *month > 2 )
    {
        *year = (unsigned)(C-4716);
    }
    else
    {
        *year = (unsigned)(C-4715);
    }
    
    return;
}

void Calendar::dateAndTimeForJulianDay(JulianDay julianDay, unsigned int *day, unsigned int *month, unsigned int *year,
                                       unsigned int *hour, unsigned int *minute, unsigned int *second)
{
    // Calculate date with float value for day.
    float floatingDay;
    Calendar::dateForJulianDay(julianDay, &floatingDay, month, year);
    
    // Store int value for day
    *day = (int)floatingDay;
    
    // Keep only decimal value
    floatingDay -= *day;
    
    // Calculate hour,minute,second
    *hour   = floatingDay*24;
    *minute = floatingDay*1440-*hour*60;
    *second = floatingDay*86400-*hour*3600-*minute*60;
    
    return;
}

unsigned int Calendar::weekDayForDate(unsigned int day, unsigned int month, unsigned int year)
{
    return weekDayForJulianDay(julianDayForDate(day,month,year));
}

unsigned int Calendar::weekDayForJulianDay(JulianDay julianDay)
{
    return (unsigned int)(julianDay.day+julianDay.time+1.5)%7;
}

