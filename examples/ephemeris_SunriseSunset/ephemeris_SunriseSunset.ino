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

void setup() 
{
  Serial.begin(9600);
  
  // Set your location on Earth
  Ephemeris::setLocationOnEarth(48,50,11,  // Lat: 48°50'11"  (Paris, France)
                                -2,20,14); // Lon: -2°20'14"

  // You can also set your altitude to improve rise and set precision
  Ephemeris::setAltitude(50);

  // Compute Sun data for a specific date and time (time is not really important for sunrise and sunset)
  SolarSystemObject sun = Ephemeris::solarSystemObjectAtDateAndTime(Sun,
                                                                    01,04,2016,
                                                                    12,00,00);

  // Print sunrise and sunset if available according to location on Earth
  if( sun.riseAndSetState == RiseAndSetOk )
  {
    int hours,minutes;
    float seconds;

    // Convert floating hours to hours, minutes, seconds and display.
    Ephemeris::floatingHoursToHoursMinutesSeconds(sun.rise, &hours, &minutes, &seconds);
    Serial.print("Sunrise: ");
    Serial.print(hours);
    Serial.print("h");
    Serial.print(minutes);
    Serial.print("m");
    Serial.print(seconds,0);
    Serial.println("s");

    // Convert floating hours to hours, minutes, seconds and display.
    Ephemeris::floatingHoursToHoursMinutesSeconds(sun.set, &hours, &minutes, &seconds);
    Serial.print("Sunset:  ");
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
}

void loop() { }
