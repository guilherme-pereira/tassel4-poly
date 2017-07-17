package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

import java.util.ArrayList;
import java.util.Vector;

/**
 * <p>Title: </p>
 * <p>Description: This will find the categories and estimate the box and wiskers for a Tablereport</p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author Ed Buckler
 * @version 1.0
 */

public class TableReportBoxWhiskerCatDataset extends DefaultBoxAndWhiskerCategoryDataset {
  String[] seriesNames;

  public TableReportBoxWhiskerCatDataset(TableReport theTable, int seriesCategory, int[] seriesY) {
    setTableReport(theTable, seriesCategory, seriesY);
  }

  public boolean setTableReport(TableReport theTable, int seriesCategory, int[] seriesY) {
    Object[][] theRawData = theTable.getTableData();
    Vector theCategories = new Vector();
    seriesNames = new String[seriesY.length];
    Object[] theSN = theTable.getTableColumnNames();
    for(int x=0; x<seriesY.length; x++){seriesNames[x]=theSN[seriesY[x]].toString();}
    for (int i = 0; i < theRawData.length; i++) {
      if(theCategories.contains(theRawData[i][seriesCategory])==false) {
        theCategories.add(theRawData[i][seriesCategory]);
      }
    }
    ArrayList[][] catData=new ArrayList[theCategories.size()][seriesY.length];
    for (int i = 0; i < theCategories.size(); i++) {
       for(int x=0; x<seriesY.length; x++){catData[i][x]=new ArrayList();}
    }
    for (int i = 0; i < theRawData.length; i++) {
      int cat=theCategories.indexOf(theRawData[i][seriesCategory]);
      Double d;
      for(int x=0; x<seriesY.length; x++){
        try {d=new Double(theRawData[i][seriesY[x]].toString());
            if(d.isNaN()==false) catData[cat][x].add(d);}
        catch (NumberFormatException ex) {}
      }
    }
    for (int i = 0; i < theCategories.size(); i++) {
      for (int x = 0; x < seriesY.length; x++) {
        //This throw errors when were are zero, but the graph is fine
        this.add(catData[i][x], seriesNames[x], theCategories.get(i).toString());
      }
    }
      return true;
  }

}
