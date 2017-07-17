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
package net.maizegenetics.util;

import java.util.Comparator;
/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */

public class StringNumberComparator implements Comparator {

  public StringNumberComparator() {
  }

  public int compare(Object o1, Object o2) {
    try{
    Double d1=new Double(o1.toString());
    Double d2=new Double(o2.toString());
    return d1.compareTo(d2);
    }
    catch(NumberFormatException e)
      {String s1=o1.toString();
      String s2=o2.toString();
      return s1.compareTo(s2);
      }

  }
  /*
  public boolean equals(Object obj) {
    throw new java.lang.UnsupportedOperationException("Method equals() not yet implemented.");
  }
*/
}