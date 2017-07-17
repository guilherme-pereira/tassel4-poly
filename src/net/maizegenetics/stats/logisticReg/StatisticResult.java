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

/**
 * Title:        TASSEL
 * Description:
 * @author Ed Buckler
 * @version 1.0
 */

public interface StatisticResult {

  public double testStatistic();

  public double testPValue();

  public double testLikelihood();

  public String[] getAttributes();

  public String getAttribute(int a);

}