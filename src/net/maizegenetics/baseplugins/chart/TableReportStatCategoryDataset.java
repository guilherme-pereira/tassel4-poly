package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;

import java.util.ArrayList;
import java.util.Vector;

/**
 * <p>Title: </p>
 * <p>Description: This will find the categories and estimate the mean and stdeviation for a Tablereport</p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author Ed Buckler
 * @version 1.0
 */

public class TableReportStatCategoryDataset extends DefaultStatisticalCategoryDataset {
  String[] seriesNames;
  boolean useStderr=false;

  public TableReportStatCategoryDataset(TableReport theTable, int seriesCategory, int[] seriesY) {
    setTableReport(theTable, seriesCategory, seriesY);
  }

  public TableReportStatCategoryDataset(TableReport theTable, int seriesCategory, int[] seriesY, boolean stderr) {
    this.useStderr=stderr;
    setTableReport(theTable, seriesCategory, seriesY);
  }


  public boolean setTableReport(TableReport theTable, int seriesCategory, int[] seriesY) {
    Object[][] theRawData = theTable.getTableData();
//    int countGood = 0;
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
            if(d.isNaN()==false) catData[cat][x].add(d);
        }
        catch (NumberFormatException ex) {}
      }
    }
    for (int i = 0; i < theCategories.size(); i++) {
      for (int x = 0; x < seriesY.length; x++) {
        Double[] d=new Double[catData[i][x].size()];
        for(int di=0; di<d.length; di++) {d[di]=(Double)catData[i][x].get(di);}
        double mean = org.jfree.data.statistics.Statistics.calculateMean(d);
        double stdev=org.jfree.data.statistics.Statistics.getStdDev(d);
        double stderr = stdev/Math.sqrt((double)d.length);
 //       System.out.println(theCategories.get(i).toString() + " = " + mean);
        if ((!Double.isNaN(mean)) && (!Double.isNaN(stdev))) {
          if(useStderr) {this.add(mean, stderr, seriesNames[x], theCategories.get(i).toString());}
          else {this.add(mean, stdev, seriesNames[x], theCategories.get(i).toString());}
        }
      }
    }
      return true;
  }

}
