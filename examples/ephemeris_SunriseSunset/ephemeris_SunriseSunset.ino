/*
 * ephemeris_SunriseSunset.ino
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
 
#include <Ephemeris.h>

void printDate(int day, int month, int year)
{
  Serial.print(day);
  Serial.print("/");
  Serial.print(month);
  Serial.print("/");
  Serial.print(year);
}

void printRiseAndSet(char *city, FLOAT latitude, FLOAT longitude, int UTCOffset, int day, int month, int year, char *ref)
{
  Ephemeris::setLocationOnEarth(latitude,longitude);
    
  Serial.print(city);
  Serial.print(" (UTC");
  if( UTCOffset >= 0 )
  {
    Serial.print("+");
  }
  Serial.print(UTCOffset);
  Serial.print(")");
  Serial.println(":");
           
  SolarSystemObject sun = Ephemeris::solarSystemObjectAtDateAndTime(Sun,
                                                                    day,month,year,
                                                                    0,0,0);

    // Print sunrise and sunset if available according to location on Earth
  if( sun.riseAndSetState == RiseAndSetOk )
  {
    int hours,minutes;
    FLOAT seconds;

    // Convert floating hours to hours, minutes, seconds and display.
    Ephemeris::floatingHoursToHoursMinutesSeconds(Ephemeris::floatingHoursWithUTCOffset(sun.rise,UTCOffset), &hours, &minutes, &seconds);
    Serial.print("  Sunrise: ");
    Serial.print(hours);
    Serial.print("h");
    Serial.print(minutes);
    Serial.print("m");
    Serial.print(seconds,0);
    Serial.println("s");

    // Convert floating hours to hours, minutes, seconds and display.
    Ephemeris::floatingHoursToHoursMinutesSeconds(Ephemeris::floatingHoursWithUTCOffset(sun.set,UTCOffset), &hours, &minutes, &seconds);
    Serial.print("  Sunset:  ");
    Serial.print(hours);
    Serial.print("h");
    Serial.print(minutes);
    Serial.print("m");
    Serial.print(seconds,0);
    Serial.println("s");
  }
  else if( sun.riseAndSetState == LocationOnEarthUnitialized )
  {
    Serial.println("You must set your location on Earth first.");
  }
  else if( sun.riseAndSetState == ObjectAlwaysInSky )
  {
    Serial.println("Sun always in sky for your location.");
  }
  else if( sun.riseAndSetState == ObjectNeverInSky )
  {
    Serial.println("Sun never in sky for your location.");
  }

  Serial.print("  Ref: ");
  Serial.println(ref);
  Serial.println();
}

void setup() 
{
  Serial.begin(9600);

  int day=10,month=2,year=2017;

  Serial.print("SunRise/SunSet at ");
  printDate(day,month,year);
  Serial.println(":\n");

  //               CITY         LATITUDE    LONGITUDE     TZ   DATE             REFERENCE
  printRiseAndSet("Paris",      48.856614,    2.3522219,  +1,  day,month,year, "sunearthtools (10/2/2017): SunRise: 08:07:00 | SunSet: 18:03:19");
  printRiseAndSet("New York",   40.7127837, -74.0059413,  -5,  day,month,year, "sunearthtools (10/2/2017): SunRise: 06:55:53 | SunSet: 17:25:09");
  printRiseAndSet("Sydney",    -33.8688197, 151.2092955, +11,  day,month,year, "sunearthtools (10/2/2017): SunRise: 06:25:25 | SunSet: 19:52:49");
  printRiseAndSet("Sao Paulo", -23.5505199, -46.6333094,  -2,  day,month,year, "sunearthtools (10/2/2017): SunRise: 06:51:35 | SunSet: 19:49:36");
}

void loop() { }
