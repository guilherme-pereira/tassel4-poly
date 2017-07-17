/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.stats.logisticReg;

/* This software is in the public domain.
   Bryan Lewis
   Dept. of Mathematics & Computer Science
   Kent State University
   url:http://www.mcs.kent.edu/~blewis
*/


public class Numeric extends Number {
/* 	This class overcomes a precision shortcoming in Java 1.02.
	All toString() conversions returned at most six
	digits. The class also provides an easy way to produce
	truncated decimal output for nice display.
*/

  public static  double MIN_VALUE = 5e-324;
  public static  double MAX_VALUE = 1.7976931348623157e+308;
  public static  double NEGATIVE_INFINITY = -1.0/0.0;
  public static  double POSITIVE_INFINITY = 1.0/0.0;
  public static  double NaN = 0.0/0.0;
  public double value;

  public Numeric(double d) {
    // constructor method
    value = d;		// double value of this object
  }

  public int intValue() {
    // conversion method to return an integer value
    return (int)value;
  }
  public long longValue() {
    // conversion method to return a long value
    return (long)value;
  }
  public float floatValue() {
    // conversion method to return a float value
    return (float)value;
  }
  public double doubleValue() {
    // conversion method to return a double value
    return value;
  }

  public String toString(int digits) {
    // this is the important method in the class
    // returns a string representation with digits digits
    int i, j;
    double t;
    char n = ' ';
    String s = new String();
    // initialize...
    t = Math.abs(value);
    j = (int) t;
    if (value < 0) {
      n = '-';
    }
    s = s + n + j + ".";
    t = t - j;
    for (i=0; i<digits; i++) {
      t = Math.abs(10 * t);
      j = (int) t;
      s = s + "" + j;
      t = t - j;
    }
    return s;
  }

}
