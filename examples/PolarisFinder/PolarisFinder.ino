/*
 * PolarisFinder.ino
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

/*
 * This sample code demonstrates how to use Ephemeris with an Arduino Mega
 * and a TFT screen to draw a Polaris finder for Takahashi mounts (P2-Z, 
 * EM-10, NJP, EM-500).
 * http://www.takahashi-europe.com/support/softwares/polarisfinder/polarisfinder-1.5-en.htm
 */
 
#include <Ephemeris.h>

#include "PolarisFinder.hpp"

#define LCD_CS A3 // Chip Select goes to Analog 3
#define LCD_CD A2 // Command/Data goes to Analog 2
#define LCD_WR A1 // LCD Write goes to Analog 1
#define LCD_RD A0 // LCD Read goes to Analog 0
#define LCD_RS A4 // LCD Reset goes to Analog 4

PolarisFinder *polarisFinder = NULL;

void setup() 
{
  Serial.begin(9600);

  float latitude  = 44.566667;
  float longitude = -6.0833333;

  int day=5,month=2,year=2017,hours=23,minutes=59,seconds=45;

  polarisFinder = new PolarisFinder(LCD_CS, LCD_CD, LCD_WR, LCD_RD, LCD_RS);
  polarisFinder->setLatitudeAndLongitude(latitude,longitude);
  polarisFinder->setDateAndTime(day,month,year,hours,minutes,seconds);
}

void loop() 
{
  polarisFinder->update();
}
