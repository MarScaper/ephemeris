# ephemeris

## Synopsis

Simple C++ library allowing to compute planet coordinates (equatorial and horizontal) with an Arduino Mega. The code is also compatible with PC.

## Features
- Compute equatorial coordinates (RA/Dec), horizontal coordinates (Alt/Az), rise, set, distance (in A.U), and apparent diameter (in arcminutes) for Sun, Mercury, Venus, Earth's Moon, Mars, Jupiter, Saturn, Uranus, Neptune.
- Uses VSOP87 for planets.
- Uses ELP2000 for the Moon.
- Code optimized for single-precision calculations with Arduino.

## Goodies
- Transformations from equatorial coordinates (RA/Dec) to horizontal coordinates (Alt/Az).
- Transformations from horizontal coordinates (Alt/Az) to equatorial coordinates (RA/Dec).
- Rise and set for custom equatorial coordinates.
- Conversion from equatorial coordinates for a specified equinox (Ex:J2000) to apparent equatorial coordinates (JNow).

## Limitations

Due to VSOP87 and ELP2000 implementation, code needs too much flash memory for classic Arduinos (Uno, etc).

## Code Example

// Set location on earth for horizontal coordinates transformations (Lat:48°50'11", Lon:-2°20'14")

Ephemeris::setLocationOnEarth(48,50,11, -2,20,14);

// Choose a date and time

int day=10,month=4,year=2014,hour=19,minute=21,second=0;

// Compute coordinates for the planet you want

SolarSystemObject planet = Ephemeris::solarSystemObjectAtDateAndTime(Mars, day, month, year, hour, minute, second);

See ephemeris_full or ephemeris_SunriseSunset in examples directory for full sample code.

## Motivation

This library is part of a personnal project to improve my EM10 Takahashi mount thanks to Arduino.

## Contributors

Code based on Jean Meeus book.

## License

Copyright (c) 2017 by Sebastien MARCHAND 
Web:www.marscaper.com - Email:sebastien@marscaper.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
