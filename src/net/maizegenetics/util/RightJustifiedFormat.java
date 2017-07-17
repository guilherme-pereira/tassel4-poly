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

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */

public class RightJustifiedFormat {

  final static double zero = 0.00001;

  public RightJustifiedFormat() {
  }

  public static String monoString(String s) {
	// convert a generic string s to a fixed length one
	String sAdd = new String(s + "                                      ");
	return sAdd.substring(0,14)+" ";
    }
    public static String monoString(double s) {
	// convert a generic double s to a "nice" fixed length string
	Double sD = new Double(s);
	String sAdd = new String();
        if(Double.isNaN(s)) {sAdd = "NaN"; }
        else if(Double.isInfinite(s)) {sAdd = "Inf";}
        else if(s>zero){sAdd = new String(sD.toString()); }
	else{ sAdd = "<0.00001"; }
	sAdd=sAdd.toLowerCase();
	int i=sAdd.indexOf('e');
	if(i>0){
	    sAdd = sAdd.substring(0,4)+"E"+sAdd.substring(i+1,sAdd.length());
	}
	else{
	    if(sAdd.length()>10){
		sAdd = sAdd.substring(0,10); }
	}
	sAdd = sAdd + "                                      ";
	return sAdd.substring(0,14)+" ";
    }
    public static String monoString(int s) {
	// convert a generic integer s to a fixed length string
	Integer sD = new Integer(s);
	String sAdd = new String(sD.toString());
	sAdd = sAdd + "                                      ";
	return sAdd.substring(0,14)+" ";
        }

    public static String monoString(Object s) {
	// convert a generic integer s to a fixed length string
        String test=s.toString();
        try{int i=Integer.parseInt(test); return monoString(i);}
        catch(Exception e) {}
        try{double d=Double.parseDouble(test); return monoString(d);}
        catch(Exception e) {}
        return monoString(test);
        }


}