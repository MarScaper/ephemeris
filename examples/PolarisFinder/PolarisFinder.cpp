/*
 * PolarisFinder.cpp
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
 
#include <Adafruit_GFX.h>    // Core graphics library
#include <Adafruit_TFTLCD.h>

#include "PolarisFinder.hpp"

#define BLACK   0x0000
#define RED     0xF800

PolarisFinder::PolarisFinder(uint8_t  LCD_CS, uint8_t LCD_CD, uint8_t LCD_WR, uint8_t LCD_RD, uint8_t LCD_RS)
{
    _tft = new Adafruit_TFTLCD(LCD_CS, LCD_CD, LCD_WR, LCD_RD, A4);
    
    _tft->reset();
    _tft->begin(_tft->readID());
    _tft->setRotation(0);
    _tft->setTextColor(RED,BLACK);
    _tft->setTextSize(2);
    _tft->setTextWrap(false);
    
    _tft->fillScreen(BLACK);

    _polarisHCoord.azi = -1000000;

    _fulldraw = true;
}

PolarisFinder::~PolarisFinder()
{
    delete _tft;
}

void PolarisFinder::setLatitudeAndLongitude(float aLatitude, float aLongitude)
{
    _latitude  = aLatitude;
    _longitude = aLongitude;

    Ephemeris::setLocationOnEarth(90,_longitude);
    
    // Flip longitude (East-/West+) to reflect PolarisFinder web app
    // provided by Optique Unterlinden for Takahashi mounts...
    Ephemeris::flipLongitude(true);
}

void PolarisFinder::setDateAndTime(int aDay, int aMonth, int aYear, int someHours, int someMinutes, int someSeconds)
{
    // Set arduino clock
    setTime(someHours, someMinutes, someSeconds, aDay, aMonth, aYear);
    
    // Set Polaris coordinates according to precession of the equinoxes
    _polarisEqCoordinates = polarisEquatorialCoordinates(aYear);
}

bool PolarisFinder::equation(float p1x, float p1y, float p2x, float p2y, float *a, float *b)
{
    if( p2x-p1x == 0 )
    {
        // [p1;p2] is vertical.
        return false;
    }
    
    *a = (p2y-p1y)/(p2x-p1x);
    *b = p1y-*a*p1x;
    
    return true;
}

EquatorialCoordinates PolarisFinder::polarisEquatorialCoordinates(long year)
{
    EquatorialCoordinates polarisEqCoord;
    
    EquatorialCoordinates coord[4];
    
    // 2000
    coord[0].ra  = 2  + 32/60.0;
    coord[0].dec = 89 + 15.9/60.0;
    
    // 2010
    coord[1].ra  = 2  + 44/60.0;
    coord[1].dec = 89 + 18.4/60.0;
    
    // 2020
    coord[2].ra  = 2  + 57/60.0;
    coord[2].dec = 89 + 20.9/60.0;
    
    // 2030
    coord[3].ra  = 3  + 12/60.0;
    coord[3].dec = 89 + 23.2/60.0;
    
    float a,b;
    
    if( year < 2010)
    {
        // 2010 and before...
        
        PolarisFinder::equation(2000,coord[0].ra,
                                2010,coord[1].ra,
                                &a,&b);
        
        polarisEqCoord.ra = a*year+b;
        
        
        PolarisFinder::equation(2000,coord[0].dec,
                                2010,coord[1].dec,
                                &a,&b);
        
        polarisEqCoord.dec = a*year+b;
    }
    else if( year >= 2010 && year < 2020 )
    {
        // Between 2010 and 2020...
        
        PolarisFinder::equation(2010,coord[1].ra,
                                2020,coord[2].ra,
                                &a,&b);
        
        polarisEqCoord.ra = a*year+b;
        
        
        PolarisFinder::equation(2010,coord[1].dec,
                                2020,coord[2].dec,
                                &a,&b);
        
        polarisEqCoord.dec = a*year+b;
    }
    else if( year >= 2020 )
    {
        // 2020 and above...
        
        PolarisFinder::equation(2020,coord[2].ra,
                                2030,coord[3].ra,
                                &a,&b);
        
        polarisEqCoord.ra = a*year+b;
        
        
        PolarisFinder::equation(2020,coord[2].dec,
                                2030,coord[3].dec,
                                &a,&b);
        
        polarisEqCoord.dec = a*year+b;
    }
    
    return polarisEqCoord;
}

void PolarisFinder::update()
{
    time_t currentTime = now();
    if( currentTime != _lastUpdate )
    {
        // One more second...

        if( _fulldraw )
        {
            drawPolarisFinder(_fulldraw);
        }
        
        int y    = 250;
        int step = 24;
        
        if( _fulldraw || day() != _day || month()!=_month || year()!=_year )
        {
            _day   = day();
            _month = month();
            _year  = year();
            
            drawDate(_day,_month,_year,60,y);
            y += step;
        }
        else
        {
            // Jump date
            y += step;
        }
        
        if( _fulldraw || hour() != _hours || minute()!=_minutes )
        {
            _hours   = hour();
            _minutes = minute();
            _seconds = second();
            
            drawTime(_hours,_minutes,70,y);
            y += step;
        }
        else
        {
            // Jump time
            y += step;
        }

        if( _fulldraw )
        {
          // Update location
          drawLatitudeOrLongitude(_latitude,0,y);
          drawLatitudeOrLongitude(_longitude,120,y);
        }

        HorizontalCoordinates newPolarisHCoord = Ephemeris::equatorialToHorizontalCoordinatesAtDateAndTime(_polarisEqCoordinates,
                                                                                                           day(), month(), year(),
                                                                                                           hour(), minute(), second());
        if( _fulldraw )
        {
            // Polaris finder needs update
            _polarisHCoord = newPolarisHCoord;
            drawPolarisAtAngle(_polarisHCoord.azi,RED);
        }
        else if( round(newPolarisHCoord.azi) != round(_polarisHCoord.azi) )
        {
            // Clear last polaris location
            drawPolarisAtAngle(_polarisHCoord.azi,BLACK);
            
            // Polaris finder needs update
            _polarisHCoord = newPolarisHCoord;
            drawPolarisFinder(_fulldraw);
            drawPolarisAtAngle(_polarisHCoord.azi,RED);
        }

        _fulldraw = false;
        
        _lastUpdate = currentTime;
    }
}

void PolarisFinder::drawLatitudeOrLongitude(float latitude, int x, int y)
{
    // Lat: 48°50'11"\0
    // Lon:  -2°20'14"\0
    char latOrLon[12];
    
    int degrees,minutes;
    float seconds;
    Ephemeris::floatingDegreesToDegreesMinutesSeconds(latitude,&degrees,&minutes,&seconds);
    
    sprintf(latOrLon,"%02d%c%02d'%02d\"",degrees,(char)247,minutes,round(seconds));
    
    _tft->setCursor(x,y);
    _tft->println(latOrLon);
}

void PolarisFinder::drawDate(int day, int month, int year, int x, int y)
{
    // 02/11/1978\0
    char date[11];
    
    sprintf(date,"%02d/%02d/%02d",day,month,year);
    
    _tft->setCursor(x,y);
    
    _tft->println(date);
}

void PolarisFinder::drawTime(int hours, int minutes, int x, int y)
{
    // 22:32 UT\0
    char time[11];
    
    sprintf(time,"%02d:%02d UT",hours,minutes);
    
    _tft->setCursor(x,y);
    
    _tft->println(time);
}

void PolarisFinder::drawPolarisFinder(bool fulldraw)
{
    int width = 12;
    
    if( fulldraw )
    {
        _tft->drawCircle(120,120,118,RED);
        _tft->drawCircle(120,120,117,RED);
        _tft->drawCircle(120,120,116,RED);
        _tft->drawCircle(120,120,116-width,RED);
    }
    
    _tft->drawCircle(120,120,116-width*2,RED);
    _tft->drawCircle(120,120,116-width*3,RED);
    
    _tft->drawLine(120,0,120,240,RED);
    _tft->drawLine(0,120,240,120,RED);
    
    for(int theta=0; theta<360; theta+=5 )
    {
        int x1 = round(80 * cos(theta*0.0174532925)+120);
        int y1 = round(80 * sin(theta*0.0174532925)+120);
        
        float r = 74;
        
        if( !fmod(theta,30) )
        {
            r = 66;
        }
        
        int x2 = round(r * cos(theta*0.0174532925)+120);
        int y2 = round(r * sin(theta*0.0174532925)+120);
        
        _tft->drawLine(x1,y1,x2,y2,RED);
    }
    
    if( fulldraw ) 
    {
        for(int theta=0; theta<360; theta+=5 )
        {
            int x1 = round(104 * cos(theta*0.0174532925)+120);
            int y1 = round(104 * sin(theta*0.0174532925)+120);
            
            float r = 108;
            
            if( !fmod(theta,30) )
            {
                r = 116;
            }
            
            int x2 = round(r * cos(theta*0.0174532925)+120);
            int y2 = round(r * sin(theta*0.0174532925)+120);
            
            _tft->drawLine(x1,y1,x2,y2,RED);
        }
    }
}

void PolarisFinder::drawPolarisAtAngle(float angle, int color)
{
    angle += 90;
    angle *= -1;
    
    int x1 = round(92 * cos(angle*0.0174532925)+120);
    int y1 = round(92 * sin(angle*0.0174532925)+120);
    
    _tft->fillCircle(x1,y1,5,color);
    
    _tft->drawLine(x1,y1,120,120,color);
}
