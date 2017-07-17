/*
 * TableReportManhattanDataset
 */
package net.maizegenetics.baseplugins.chart;

import java.util.HashMap;
import net.maizegenetics.pal.report.TableReport;
import org.jfree.data.xy.DefaultTableXYDataset;

/**
 *
 * @author yz79
 */
public class TableReportManhattanDataset extends DefaultTableXYDataset {

    double[][] theData;
    String[] seriesNames;
    String xName;
    String myTrait;
    int numberYAxes;
    String[] myChromNames;
    Object[] myColumnNames;
    double[] myPValues;
    double[] myLogPValues;
    int[] myPositions;
    String[] myMarkers;
    HashMap myLookupTable;
    int myPValueColumnIndex = -1;
    int myChromColumnIndex = -1;
    int myPositionColumnIndex = -1;
    int myMarkerColumnIndex = -1;
    int myTraitColumnIndex = -1;
    int myNumRows;
    int myStartIndex;
    int myEndIndex;
    boolean myNumericChromNames = true;

    public TableReportManhattanDataset(TableReport theTable, int start, int end) {
        numberYAxes = 0;
        myStartIndex = start;
        myEndIndex = end;
        myNumRows = myEndIndex - start;
        setTableReport(theTable);
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

    private void setTraitColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Trait")) {
                myTraitColumnIndex = i;
                return;
            }
        }
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

    private void setChromColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Locus")) {
                myChromColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No Chromosome names in selected data");
    }

    private void setMarkerColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Marker")) {
                myMarkerColumnIndex = i;
                return;
            }
        }
    }

    private void setPositionColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Locus_pos") || myColumnNames[i].equals("Site")) {
                myPositionColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No positions in selected data");
    }

    private void setChromNames(TableReport myTableReport) {
        String currentChrom = "";
        for (int i = 0; i < myChromNames.length; i++) {
            myChromNames[i] = ((String) myTableReport.getValueAt(myStartIndex + i, myChromColumnIndex));
            if (!currentChrom.equals(myChromNames[i])) {
                numberYAxes++;
                currentChrom = myChromNames[i];
            }
        }
    }

    private void setNumericChromNames() {
        for (int i = 0; i < myChromNames.length; i++) {
            try {
                if (numberYAxes < Integer.parseInt(myChromNames[i])) {
                    numberYAxes = Integer.parseInt(myChromNames[i]);
                }
            } catch (Exception e) {
                myNumericChromNames = false;
                break;
            }
        }
    }

    private void setMarkers(TableReport myTableReport) {
        for (int i = 0; i < myMarkers.length; i++) {
            myMarkers[i] = ((String) myTableReport.getValueAt(myStartIndex + i, myMarkerColumnIndex));
        }
    }

    private void setPValues(TableReport myTableReport) {
        Object temp = myTableReport.getValueAt(myStartIndex, myPValueColumnIndex);
        if (temp instanceof Double) {
            for (int i = 0; i < myPValues.length; i++) {
                myPValues[i] = ((Double) myTableReport.getValueAt(myStartIndex + i, myPValueColumnIndex)).doubleValue();
                if (myPValues[i] == 0) {
                    myPValues[i] = Double.MIN_VALUE;
                }
                myLookupTable.put(myPValues[i], i);
            }
        } else if (temp instanceof String) {
            for (int i = 0; i < myPValues.length; i++) {
                myPValues[i] = Double.parseDouble((String) myTableReport.getValueAt(myStartIndex + i, myPValueColumnIndex));
                if (myPValues[i] == 0) {
                    myPValues[i] = Double.MIN_VALUE;
                }
                myLookupTable.put(myPValues[i], i);
            }
        } else {
            throw new IllegalStateException("TableReportManhattanDataset: setPValues: Unknown data type of P values: " + temp.getClass().getName());
        }
    }

    private void setPositions(TableReport myTableReport) {
        int offset = 0;
        int previousPosition = 0;
        int currentPosition = 0;
        // GLM positions formated as int
        for (int i = 0; i < myPositions.length; i++) {
            currentPosition = Integer.valueOf((myTableReport.getValueAt(myStartIndex + i, myPositionColumnIndex)).toString());
            myPositions[i] = currentPosition + offset;
            if (currentPosition < previousPosition) {
                offset = offset + previousPosition;
            }
            previousPosition = currentPosition;
        }
    }

    private void setTrait(TableReport table) {
        myTrait = (String) table.getValueAt(myStartIndex, myTraitColumnIndex);
    }

    public String[] getMarkers() {
        return myMarkers;
    }

    public String[] getChroms() {
        return myChromNames;
    }

    public String getTrait() {
        return myTrait;
    }

    private void setLogPValues() {
        for (int i = 0; i < myLogPValues.length; i++) {
            myLogPValues[i] = -Math.log10(myPValues[i]);
        }
    }

    public void setTableReport(TableReport theTable) {
        myColumnNames = theTable.getTableColumnNames();
        setPValueColumnIndex();
        setChromColumnIndex();
        setPositionColumnIndex();
        setMarkerColumnIndex();
        setTraitColumnIndex();
        myPValues = new double[myNumRows];
        myLogPValues = new double[myNumRows];
        myChromNames = new String[myNumRows];
        myPositions = new int[myNumRows];
        myMarkers = new String[myNumRows];
        myLookupTable = new HashMap(myNumRows);
        setPValues(theTable);
        setLogPValues();
        setChromNames(theTable);
        setNumericChromNames();
        setPositions(theTable);
        setMarkers(theTable);
        setTrait(theTable);
        seriesNames = new String[numberYAxes];
        seriesNames[0] = myChromNames[0];

        String currentChrom = myChromNames[0];
        int chromIndex = 1;
        theData = new double[myNumRows][1 + numberYAxes];
        for (int i = 0; i < theData.length; i++) {
            for (int j = 0; j < theData[0].length; j++) {
                theData[i][j] = Double.NaN;
            }
        }
        for (int i = 0; i < myNumRows; i++) {
            try {
                theData[i][0] = myPositions[i];
                if (!myNumericChromNames) {
                    if (!currentChrom.equals(myChromNames[i])) {
                        chromIndex++;
                        currentChrom = myChromNames[i];
                        seriesNames[chromIndex - 1] = currentChrom;
                    }
                    //                ÷if (!myNumericChromNames) {
                    theData[i][chromIndex] = myLogPValues[i];
                } else {

                    for (int j = 0; j < numberYAxes; j++) {
                        seriesNames[j] = Integer.toString(j + 1);
                    }

                    theData[i][Integer.parseInt(myChromNames[i])] = myLogPValues[i];

                }
//               ÷ } else {
//                    theData[i][myChromOrder[chromIndex - 1]+1] = myLogPValues[i];
//                }
            } catch (NumberFormatException ex) {
                System.out.println("throw new NumberFormatException();");
            }
        }
        xName = "Position";
    }
}