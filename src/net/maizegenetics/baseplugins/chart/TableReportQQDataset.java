/*
 * TableReportQQDataset
 */
package net.maizegenetics.baseplugins.chart;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import net.maizegenetics.pal.report.TableReport;
import org.jfree.data.xy.DefaultTableXYDataset;

/**
 *
 * @author yz79
 */
public class TableReportQQDataset extends DefaultTableXYDataset {

    double[][] theData;
    String[] seriesNames;
    String xName;
    String myTrait;
    int numberYAxes;
    Object[] myColumnNames;
    double[] myPValues;
    double[] myLogPValues;
    double[] myExpectedPValues;
    double[] myLogExpectedPValues;
    int[] myPositions;
    String[] myMarkers;
    HashMap myLookupTable;
    int myPValueColumnIndex = -1;
    int myPositionColumnIndex = -1;
    int myTraitColumnIndex = -1;
    int myMarkerColumnIndex = -1;
    int myNumRows;
    int myStartIndex;
    int myEndIndex;
    int myCountToDisplay;
    double myDistance = 0.01;

    public TableReportQQDataset(TableReport table) {
        numberYAxes = 1;
        setTableReport(table);
    }

    public TableReportQQDataset(TableReport table, int startIndex, int endIndex, int countToDisplay) {
        numberYAxes = 1;
        myStartIndex = startIndex;
        myEndIndex = endIndex;
        myNumRows = myEndIndex - startIndex;
        myCountToDisplay = countToDisplay;
//        myNumRows = countToDisplay + (myNumRows - countToDisplay)/20 + 1;
        setTableReport(table);
    }

    public int getItemCount(int parm1) {
        return theData.length;
        //throw new java.lang.UnsupportedOperationException("Method getItemCount() not yet implemented.");
    }

    public Number getX(int series, int item) {
        Double x = new Double(theData[item][0]);
        return x;
        //    throw new java.lang.UnsupportedOperationException("Method getXValue() not yet implemented.");
    }

    public int getSeriesCount() {
        return numberYAxes;
        //throw new java.lang.UnsupportedOperationException("Method getSeriesCount() not yet implemented.");
    }

    public Number getY(int series, int item) {
        Double y = new Double(theData[item][1 + series]);
        return y;
        //    throw new java.lang.UnsupportedOperationException("Method getYValue() not yet implemented.");
    }

    public String getSeriesName(int series) {
        /**
         * current
         */
        return seriesNames[series];
        //    throw new java.lang.UnsupportedOperationException("Method getSeriesName() not yet implemented.");
    }

    public String getSeriesKey(int series) {
        /**
         * current
         */
        return seriesNames[series];
        //    throw new java.lang.UnsupportedOperationException("Method getSeriesName() not yet implemented.");
    }

    public String getXName() {
        /**
         * current
         */
        return xName;
        //    throw new java.lang.UnsupportedOperationException("Method getSeriesName() not yet implemented.");
    }

    private void setPValueColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("p") || myColumnNames[i].equals("marker_p")) {
                myPValueColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No P-values in selected data");
    }

    private void setPositionColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Locus_pos") || myColumnNames[i].equals("Site")) {
                myPositionColumnIndex = i;
                return;
            }
        }
    }

    private void setTraitColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Trait")) {
                myTraitColumnIndex = i;
                return;
            }
        }
    }

    private void setMarkerColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Marker")) {
                myMarkerColumnIndex = i;
                return;
            }
        }
    }

    private void setTrait(TableReport table) {
        myTrait = (String) table.getValueAt(myStartIndex, myTraitColumnIndex);
    }

    private void setPValues(TableReport myTableReport) {
        Object temp = myTableReport.getValueAt(myStartIndex, myPValueColumnIndex);
        if (temp instanceof Double) {
            for (int i = 0; i < myPValues.length; i++) {
                myPValues[i] = ((Double) myTableReport.getValueAt(myStartIndex + i, myPValueColumnIndex)).doubleValue();
                if (myPValues[i] == 0) {
                    myPValues[i] = Double.MIN_VALUE;
                }
                myLookupTable.put(-Math.log10(myPValues[i]), i);
            }
        } else if (temp instanceof String) {
            for (int i = 0; i < myPValues.length; i++) {
                myPValues[i] = Double.parseDouble((String) myTableReport.getValueAt(myStartIndex + i, myPValueColumnIndex));
                if (myPValues[i] == 0) {
                    myPValues[i] = Double.MIN_VALUE;
                }
                myLookupTable.put(-Math.log10(myPValues[i]), i);
            }
        } else {
            throw new IllegalStateException("TableReportQQDataset: setPValues: Unknown data type of P values: " + temp.getClass().getName());
        }
    }

    private void setPositions(TableReport myTableReport) {
//        if (myColumnNames[myPositionColumnIndex].equals("Locus_pos")) {
//            for (int i = 0; i < myPositions.length; i++) {
//                myPositions[i] = ((Integer)myTableReport.getValueAt(myStartIndex + i, myPositionColumnIndex)).intValue();
//            }
//        }
//        else if (myColumnNames[myPositionColumnIndex].equals("Site")) {
        for (int i = 0; i < myPositions.length; i++) {
            myPositions[i] = Integer.valueOf((myTableReport.getValueAt(myStartIndex + i, myPositionColumnIndex)).toString());
        }
//        }
    }

    private void setMarkers(TableReport myTableReport) {
        for (int i = 0; i < myMarkers.length; i++) {
            myMarkers[i] = ((String) myTableReport.getValueAt(myStartIndex + i, myMarkerColumnIndex));
        }
    }

    public int[] getPositions() {
        return myPositions;
    }

    public Map getLookupTable() {
        return myLookupTable;
    }

    public String[] getMarkers() {
        return myMarkers;
    }

    private void sortPValues() {
        Arrays.sort(myPValues);
    }

    private void setLogPValues() {
        for (int i = 0; i < myLogPValues.length; i++) {
            myLogPValues[i] = -Math.log10(myPValues[i]);
        }
    }

    private void setExpectedPValues() {
        double increment = 1 / (double) myNumRows;
        for (int i = 0; i < myExpectedPValues.length; i++) {
            myExpectedPValues[i] = increment * ((double) (i + 1));
        }
    }

    private void setLogExpectedPValues() {
        for (int i = 0; i < myLogExpectedPValues.length; i++) {
            myLogExpectedPValues[i] = -Math.log10(myExpectedPValues[i]);
        }
    }

    public void setTableReport(TableReport theTable) {
        myColumnNames = theTable.getTableColumnNames();
        setPValueColumnIndex();
        setPositionColumnIndex();
        setMarkerColumnIndex();
        myPValues = new double[myNumRows];
        myLogPValues = new double[myNumRows];
        myExpectedPValues = new double[myNumRows];
        myLogExpectedPValues = new double[myNumRows];
        myPositions = new int[myNumRows];
        myMarkers = new String[myNumRows];
        myLookupTable = new HashMap(myNumRows);
        setPValues(theTable);
        setPositions(theTable);
        setMarkers(theTable);
        sortPValues();
        setLogPValues();
        setExpectedPValues();
        setLogExpectedPValues();
        setTraitColumnIndex();
        setTrait(theTable);
        theData = new double[myNumRows][2];
        double previousXValue = myLogExpectedPValues[0];
        for (int i = 0; i < myNumRows; i++) {
            try {
                theData[i][0] = myLogExpectedPValues[i];
                if (i < myCountToDisplay) {
                    if (myLogPValues[i] != 0) {
                        theData[i][1] = myLogPValues[i];
                    } else {
                        theData[i][1] = Double.NaN;
                    }
                } else {
                    double value = myLogExpectedPValues[i];
                    if (previousXValue - value > myDistance) {
                        theData[i][1] = myLogPValues[i];
                        previousXValue = myLogExpectedPValues[i];
                    } else {
                        theData[i][1] = Double.NaN;
                    }
                }
            } catch (NumberFormatException ex) {
                System.out.println("throw new NumberFormatException();");
            }
        }

        theData[0][0] = myLogExpectedPValues[0];
        theData[myNumRows - 1][0] = myLogExpectedPValues[myNumRows - 1];
        seriesNames = new String[1];
        xName = "Expected -Log(P-Value)";
        seriesNames[0] = myTrait;
    }
}