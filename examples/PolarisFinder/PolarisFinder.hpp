/*
 * PolarisFinder.hpp
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

#ifndef PolarisFinder_h
#define PolarisFinder_h

#include <Arduino.h>
#include <Ephemeris.h>
#include <TimeLib.h>

class Adafruit_TFTLCD;

class PolarisFinder
{
public:
    
    PolarisFinder(uint8_t  LCD_CS, uint8_t LCD_CD, uint8_t LCD_WR, uint8_t LCD_RD, uint8_t LCD_RS);
    
    ~PolarisFinder();
    
    void setLatitudeAndLongitude(float aLatitude, float aLongitude);
    
    void setDateAndTime(int aDay, int aMonth, int aYear, int someHours, int someMinutes, int someSeconds);
    
    void update();
    
private:

    bool _fulldraw;
    
    float _latitude;
    float _longitude;
    
    int _day;
    int _month;
    int _year;
    
    int _hours;
    int _minutes;
    int _seconds;
    
    // Last update time in seconds
    time_t _lastUpdate;
    
    // Screen for display
    Adafruit_TFTLCD *_tft;
    
    // Current Polaris coordinates
    EquatorialCoordinates _polarisEqCoordinates;
    HorizontalCoordinates _polarisHCoord;
    
    // Polaris coordinates according to precession of the Equinoxes
    static EquatorialCoordinates polarisEquatorialCoordinates(long year);
    
    // Used for linear interpolation of polaris coordinates
    static bool equation(float p1x, float p1y, float p2x, float p2y, float *a, float *b);
    
    // Draw coordinates on Earth
    void drawLatitudeOrLongitude(float latitude, int x, int y);
    
    // Draw date
    void drawDate(int day, int month, int year, int x, int y);
    
    // Draw time
    void drawTime(int hours, int minutes, int x, int y);
    
    // Draw Polaris finder circles and graduation
    void drawPolarisFinder(bool fulldraw);
    
    // Draw Polaris
    void drawPolarisAtAngle(float angle, int color);
};

#endif
