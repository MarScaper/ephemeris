/*
 * Calendar.hpp
 */

#ifndef Calendar_h
#define Calendar_h

/*! This structure describes a Julian day with integer and floating value parts.
 *  As a result we can get 64bits precision with single precision processor. */
struct JulianDay
{
    /*! Integer value of Jullian day. */
    long day;
    
    /*! Floating value of Jullian day (0.0=0h 1.0=24h) */
    float time;
};


/*!
 * This class manipulates calendar arithmetic.
 */
class Calendar
{
public:
    
    static JulianDay julianDayForDate(float day, unsigned int month, unsigned int year);
    
    static JulianDay julianDayForDateAndTime(unsigned int day, unsigned int month, unsigned int year,
                                             unsigned int hour, unsigned int minute, unsigned int second);
    
    static void dateForJulianDay(JulianDay julianDay, float *day, unsigned int *month, unsigned int *year);
    
    static void dateAndTimeForJulianDay(JulianDay julianDay,
                                        unsigned int *day,  unsigned int *month,  unsigned int *year,
                                        unsigned int *hour, unsigned int *minute, unsigned int *second);
    
    static unsigned int weekDayForDate(unsigned int day, unsigned int month, unsigned int year);
    
    static unsigned int weekDayForJulianDay(JulianDay julianDay);
};

#endif

