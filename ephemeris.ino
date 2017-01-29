/*
 * ephemeris.ino
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
    
#include "Ephemeris.hpp"

void printDateAndTime(int day, int month, int year, int hour, int minute, int second )
{
  Serial.print(day);
  Serial.print("/");
  Serial.print(month);
  Serial.print("/");
  Serial.print(year);
  Serial.print(" ");
  Serial.print(hour);
  Serial.print(":");
  Serial.print(minute);
  Serial.print(":");
  Serial.print(second);
}

void equatorialCoordinatesToString(EquatorialCoordinates coord, char raCoord[14] , char decCoord[14])
{
  int raHour,raMinute;
  float raSecond;
  Ephemeris::floatingHoursToHoursMinutesSeconds(coord.ra, &raHour, &raMinute, &raSecond);
    
  sprintf(raCoord," %02dh%02dm%02ds.%02d",raHour,raMinute,(int)raSecond,(int)round(((float)(raSecond-(int)raSecond)*pow(10,2))));
    
  int decDegree,decMinute;
  float decSecond;
  Ephemeris::floatingDegreesToDegreesMinutesSeconds(coord.dec, &decDegree, &decMinute, &decSecond);
    
  if(decDegree<0)
  {
    sprintf(decCoord,"%02dd%02d'%02d\".%02d",(int)decDegree,decMinute,(int)decSecond,(int)round(((float)(decSecond-(int)decSecond)*pow(10,2))));
  }
  else
  {
    sprintf(decCoord," %02dd%02d'%02d\".%02d",(int)decDegree,decMinute,(int)decSecond,(int)round(((float)(decSecond-(int)decSecond)*pow(10,2))));
  }
}

void printEquatorialCoordinates(EquatorialCoordinates coord)
{
  if( isnan(coord.ra) ||  isnan(coord.dec))
  {
    // Do not work for Earth of course...
    Serial.println("R.A: -");
    Serial.println("Dec: -");
        
    return;
  }
    
  char raCoord[14];
  char decCoord[14];
  equatorialCoordinatesToString(coord,raCoord,decCoord);

  Serial.print("R.A: ");
  Serial.println(raCoord);

  Serial.print("Dec: ");
  Serial.println(decCoord);

  return;
}

void printHorizontalCoordinates(HorizontalCoordinates coord)
{
  if( isnan(coord.azi) ||  isnan(coord.alt))
  {
    // Do not work for Earth of course...
    Serial.println("Azi: -");
    Serial.println("Alt: -");
        
    return;
  }

  Serial.print("Azi: ");
  Serial.print(coord.azi,2);
  Serial.println("d");

  Serial.print("Alt: ");
  Serial.print(coord.alt,2);
  Serial.println("d");
}

void printSolarSystemObjects(int day, int month, int year, int hour, int minute, int second)
{
  Serial.println("_____________________________________");
  printPlanet("Sun",     Sun,     day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Mercury", Mercury, day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Venus",   Venus,   day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Earth",   Earth,   day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Mars",    Mars,    day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Jupiter", Jupiter, day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Saturn",  Saturn,  day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Uranus",  Uranus,  day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
  printPlanet("Neptune", Neptune, day, month, year, hour, minute, second);
  Serial.println("_____________________________________");
}

void printPlanet(char *solarSystemObjectName, SolarSystemObjectIndex index, int day, int month, int year, int hour, int minute, int second )
{
  SolarSystemObject planet = Ephemeris::solarSystemObjectAtDateAndTime(index, day, month, year, hour, minute, second);

  Serial.println(solarSystemObjectName);
  printEquatorialCoordinates(planet.equaCoordinates);
  printHorizontalCoordinates(planet.horiCoordinates);

  if( isnan(planet.diameter) )
  {
    // Do not work for Earth of course...
    Serial.println("Dist: -");
    Serial.println("Diam: -");
  }
  else
  {
    Serial.print("Dist: ");
    Serial.print(planet.distance,3);
    Serial.println(" AU");
  
    if( planet.diameter <= 1 )
    {
      Serial.print("Diam: ");
      Serial.print(planet.diameter*60,2);
      Serial.println("\"");
    }
    else
    {
      Serial.print("Diam: ");
      Serial.print(planet.diameter,2);
      Serial.println("'");
    }
  }
}

void setup() 
{
  Serial.begin(9600);

  // Set location on earth for horizontal coordinates transformations
  Ephemeris::setLocationOnEarth(48,50,11,  // Lat: 48°50'11"
                                -2,20,14); // Lon: -2°20'14"
  // Choose a date and time
  int day=10,month=4,year=2014,hour=19,minute=21,second=0;

  // Compute and print solar system objects
  Serial.print("Coordinates of Solar system objects (");
  printDateAndTime(day,month,year,hour,minute,second);
  Serial.println(")");
  printSolarSystemObjects(day, month, year, hour, minute, second);
  
  return;
}

void loop() {}
