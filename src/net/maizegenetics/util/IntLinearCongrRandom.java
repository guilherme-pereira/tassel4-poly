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

//Title:        QTPAnalyzer
//Version:
//Copyright:    Copyright (c) 1997
//Author:       Ed Buckler
//Company:      NCSU
//Description:  Your description


package net.maizegenetics.util;
import java.util.Date;
public class IntLinearCongrRandom { //This series is fast and even but unknown # of repeats

	private int x0;
	public IntLinearCongrRandom()
		{
			long time;

			Date s=new Date();
			time=s.getTime();
			x0=(int)(time%894407);
		}

	public IntLinearCongrRandom(int seed)
		{
			x0=seed;
		}

  public void setSeed(int seed)
		{
			x0=seed;
		}

	public int nextInt(int range)
		{
			x0=((2401*x0)%894407);
			return (x0%range);
		}
}