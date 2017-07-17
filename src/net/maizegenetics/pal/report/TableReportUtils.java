package net.maizegenetics.pal.report;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

import net.maizegenetics.util.DoubleFormat;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

/**
 * @author terry
 */
public class TableReportUtils {

    public static void saveDelimitedTableReport(TableReport theTableSource, String delimit, File saveFile) {

        if (saveFile == null) {
            return;
        }

        BufferedWriter bw = null;
        try {

            bw = Utils.getBufferedWriter(saveFile);

            Object[] colNames = theTableSource.getTableColumnNames();
            for (int j = 0; j < colNames.length; j++) {
                if (j != 0) {
                    bw.write(delimit);
                }
                bw.write(colNames[j].toString());
            }
            bw.write("\n");

            for (int r = 0, n = theTableSource.getRowCount(); r < n; r++) {
                Object[] theRow = theTableSource.getRow(r);
                for (int i = 0; i < theRow.length; i++) {
                    if (i != 0) {
                        bw.write(delimit);
                    }
                    if (theRow[i] == null) {
                        // do nothing
                    } else if (theRow[i] instanceof Double) {
                        bw.write(DoubleFormat.format((Double) theRow[i]));
                    } else {
                        bw.write(theRow[i].toString());
                    }
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("TableReportUtils: writeReport: problem writing file: " + e.getMessage());
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    public static TableReport readDelimitedTableReport(String saveFile, String delimit) {

        int numLines = -1;
        BufferedReader reader = null;
        try {
            reader = Utils.getBufferedReader(saveFile, 1000000);
            while (reader.readLine() != null) {
                numLines++;
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating TableReport: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                reader.close();
            } catch (Exception ex) {
                // do nothing
            }
        }

        BufferedReader br = null;
        try {
            br = Utils.getBufferedReader(saveFile);
            String[] columnHeaders = br.readLine().trim().split(delimit);

            String[][] data = new String[numLines][];
            for (int i = 0; i < numLines; i++) {
                data[i] = br.readLine().trim().split(delimit);
            }
            return new SimpleTableReport(saveFile, columnHeaders, data);
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating TableReport: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                br.close();
            } catch (Exception ex) {
                // do nothing
            }
        }

    }
}
