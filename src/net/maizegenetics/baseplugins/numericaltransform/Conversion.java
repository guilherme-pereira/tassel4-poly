package net.maizegenetics.baseplugins.numericaltransform;

import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.alignment.Trait;

import java.math.BigDecimal;
import java.util.ArrayList;

/**
 * User: dkroon
 * Date: Aug 1, 2005
 * Time: 10:41:05 AM
 */
public class Conversion {

    /**
     *
     * @param tableReport
     * @param colSelected Which column should be parsed and returned as a double[].
     * @return
     */
    public static double[] parseColumnData(Phenotype tableReport, int colSelected) {

        double[][] rawData = tableReport.getData();
        //Object[][] rawData = tableReport.getTableData();
        double[] tempData = new double[rawData.length];
        int naNCount = 0;  // how many NaN are in the column
        for (int i = 0; i < rawData.length; i++) {     // start with 1 because taxa names are in column 0
            //tempData[i] = Double.valueOf(rawData[i][colSelected].toString()).doubleValue();
            tempData[i] = rawData[i][colSelected];
            if (Double.isNaN(tempData[i])) {
                naNCount++;
            }
        }

        // if there is not a single numeric value in the column, then go no further
        if (naNCount == rawData.length) {
            return null;
        }

        return tempData;
    }

    /**
     * Does a simple conversion of the values in a TableReport into a double[][].
     * @param theCharacterAlignment
     * @param selectedCol Columns from the TableReport which are to be included in the returned double[][].  If
     *          selectedCols == null, then it includes all columns.
     * @return double[][] of the values in the TableReport
     */
    public static double[][] parseColumnData(Phenotype theCharacterAlignment, int[] selectedCol) {
        //todo change to character alignment as tablereports aren't necessarily numeric
        // create the appropriate array if the selectedCol == null
        int[] includedColumn;
        if (selectedCol == null) {
            int colCount = theCharacterAlignment.getNumberOfTraits();
            includedColumn = new int[colCount];
            for (int i = 0; i < colCount; i++) {
                includedColumn[i] = i;
            }
        } else {
            includedColumn = selectedCol;
        }
        // do the simple conversion
        double[][] tempData = new double[theCharacterAlignment.getNumberOfTaxa()][includedColumn.length];
        for (int i = 0; i < theCharacterAlignment.getNumberOfTaxa(); i++) {
            for (int j = 0; j < includedColumn.length; j++) {
                tempData[i][j] = theCharacterAlignment.getData(i, includedColumn[j]);
            }
        }
        return tempData;
    }

    /**
     * 2D array version for iterating through the full raw data set and replaces the
     * original data with the changed data for the appropriate columns.
     * 
     * @param originalSca
     * @param colSelected
     * @param changedData
     * @return
     */
    public static SimplePhenotype reconstituteDataset(Phenotype originalSca, int[] colSelected, double[][] changedData) {
    	
    	ArrayList<Trait> newtraits = new ArrayList<Trait>();
        for (int col : colSelected) {
            newtraits.add(Trait.getInstance(originalSca.getTrait(col)));
        }
        
        return new SimplePhenotype(originalSca.getTaxa(), newtraits, changedData);
    }

    /**
     * Remove NaN values and resize the array.
     *
     * @param data Array containing NaN values.
     * @return Resized array without any NaN values.
     */
    public static double[] removeNaNValues(double[] data) {
        if (data == null || data.length == 0) {
            return null;
        }

        double[] cleanData = new double[data.length];
        int count = 0;
        for (int i = 0; i < data.length; i++) {
            if (!Double.isNaN(data[i])) {
                cleanData[count++] = data[i];
            }
        }
        double[] tempData = new double[count];
        System.arraycopy(cleanData, 0, tempData, 0, count);
        return tempData;
    }

    /**
     * Remove NaN values and resize the array.
     *
     * @param data 2-D array containing NaN values.
     * @return Resized array without any NaN values.
     */
    public static double[][] removeNaNValues(double[][] data) {
        if (data == null || data.length == 0) {
            return null;
        }
        double[][] cleanData = new double[data.length][];
        int count = 0;
        for (int i = 0; i < data.length; i++) {
            cleanData[i] = removeNaNValues(data[i]);
        }

        return cleanData;
    }

    /**
     * Determine the percentage of data which is NaN in the passed-in TableReport
     * @param tableReportIn
     * @param precision The number of desired significan digits after the decimal,
     *          i.e., the significand or, more informally, the mantissa
     * @return
     */
    public static Object[] getPercentMissingData(Phenotype tableReportIn, int precision) {

        double[][] rawData = tableReportIn.getData();
        double[] tempData = new double[rawData.length];


        int colCount = rawData[0].length;   // assume that the first row has all columns
        int rowCount = rawData.length;
        BigDecimal hundred = new BigDecimal("100");
        BigDecimal[] percentData = new BigDecimal[colCount];
        for (int j = 0; j < colCount; j++) {
            int naNCount = 0;  // how many NaN are in the column
            for (int i = 0; i < rowCount; i++) {
                tempData[i] = Double.valueOf(String.valueOf(rawData[i][j])).doubleValue();
                if (Double.isNaN(tempData[i])) {
                    naNCount++;
                }
            }
            BigDecimal bd = new BigDecimal((double) naNCount / rowCount);
            BigDecimal percentage = bd.multiply(hundred);
            percentData[j] = percentage.setScale(precision, BigDecimal.ROUND_HALF_UP);
        }

        return percentData;
    }

    /**
     * Tests the normality of the given data and returns the pre-formatted results
     * from Kolmogorov-Smirnov, Cramer-vonMises and Anderson-Darling.
     * Utilizes the Stochastic Simulation in Java (SSJ) library.
     *
     * @param data the data set to be tested
     * @return pre-formatted String of the normality test results
     */
    /*   public static String getNormalityResults(double[] data) {
    DoubleArrayList dal = new DoubleArrayList(data);

    ContinuousEmpiricalDist ced = new ContinuousEmpiricalDist(data);

    GofFormat.activeTests[2] = true;      // set Kolmogorov-Smirnov to be done
    GofFormat.activeTests[4] = true;      // set Cramer-vonMises to be done

    double[] results = new double[GofFormat.NTESTTYPES];
    double[] pResults = new double[GofFormat.NTESTTYPES];
    GofFormat.activeTests(dal, ced, results, pResults);

    return GofFormat.formatActiveTests(dal.size(), results, pResults);
    }
     */
    public static double[] normalizeData(double[] dataIn) {

        // translate into a 1-based array to prevent division by zero
        double[] data = new double[dataIn.length + 1];
        System.arraycopy(dataIn, 0, data, 1, dataIn.length);
        double avg = 0;
        double cumulativeValue = 0;
        int n = 0;
        for (int i = 1; i < data.length; i++) {
            if (!Double.isNaN(data[i])) {
                cumulativeValue += data[i];
                n++;
            }
        }
        avg = cumulativeValue / n;
        double stDev = calculateStandardDeviation(dataIn);
        double[] result = new double[dataIn.length];
        for (int i = 0; i < dataIn.length; i++) {
            if (!Double.isNaN(dataIn[i])) {
                result[i] = (dataIn[i] - avg) / stDev;
            } else {
                result[i] = Double.NaN;
            }
        }
        return result;
    }

    /**
     * Normalizes the data in each column separately.
     *
     * @param dataIn
     * @return
     */
    public static double[][] normalizeData(double[][] dataIn) {
        double[][] result = new double[dataIn.length][dataIn[0].length];

        int colCount = dataIn[0].length;
        // do all of the calculations a column at a time
        for (int j = 0; j < colCount; j++) {
            double[] colData = new double[dataIn.length];
            for (int row = 0; row < dataIn.length; row++) {
                colData[row] = dataIn[row][j];
            }
            double[] normalizedData = normalizeData(colData);

            // fill the rows into the same column
            for (int q = 0; q < dataIn.length; q++) {
                if (!Double.isNaN(dataIn[q][j])) {
                    result[q][j] = normalizedData[q];
                } else {
                    result[q][j] = Double.NaN;
                }
            }
        }
        return result;
    }

    public static double calculateMean(double[] data) {
        if (data.length == 0) {
            return 0;
        }
        double sum = 0;
        for (int i = 0; i < data.length; i++) {
            if (!Double.isNaN(data[i])) {
                sum += data[i];
            }
        }
        return (double) sum / data.length;
    }

    /**
     * Calculate variance for a single array of values.
     *
     * @param dataIn
     * @return
     */
    public static double calculateVariance(double[] dataIn) {
        // translate into a 1-based array to prevent division by zero
        double[] data = new double[dataIn.length + 1];
        System.arraycopy(dataIn, 0, data, 1, dataIn.length);

        int n = 0;
        double sum = 0;
        double sumSqr = 0;

        for (int i = 1; i < data.length; i++) {

            if (!Double.isNaN(data[i])) {
                n += 1;
                sum += data[i];
                sumSqr += data[i] * data[i];
            }
        }

        return (sumSqr - sum * sum / n) / (n - 1);
    }

    /**
     *
     * @param data
     * @return
     */
    public static double calculateStandardDeviation(double[] data) {
        return Math.sqrt(calculateVariance(data));
    }
}
