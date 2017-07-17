package net.maizegenetics.stats.GLM;

//import net.maizegenetics.stats.Statistic;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.pal.report.AbstractTableReport;
import net.maizegenetics.pal.report.TableReport;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Apr 26, 2004
 * Time: 12:15:26 PM
 * <p/>
 * Report writer handles Statistics with one or two names in the name array.
 * It calls the first name factor and the second level
 */
public class ReportWriter extends AbstractTableReport implements TableReport, Serializable {
    //fields
//    private Object[][] contents;

    ArrayList<Object[]> contents;
    private String[] colNames;
    private String title;

    //constructors
    public ReportWriter(String[] columnNames, Object[][] results, String title) {

        colNames = columnNames;
        contents = new ArrayList<Object[]>();
        for (int row = 0; row < results.length; row++) {
            contents.add(results[row]);
        }
        this.title = title;
        ;
    }

    public ReportWriter(String[] columnNames, String title) {
        colNames = columnNames;
        contents = new ArrayList<Object[]>();
        this.title = title;
    }

    public ReportWriter(Reporter rp, String title) {
        colNames = rp.getColumnNames();
        Object[][] results = rp.getResults();
        contents = new ArrayList<Object[]>();
        for (int row = 0; row < results.length; row++) {
            contents.add(results[row]);
        }

        this.title = title;

    }


    //methods
    public void append(Object[][] results) {
        for (int i = 0; i < results.length; i++) {
            contents.add(results[i]);
        }
    }

    public void append(Reporter rp) {
        append(rp.getResults());
    }

    public Object[] getTableColumnNames() {
        return colNames;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {

        DecimalFormat dfexp = new DecimalFormat("0.####E0");
        dfexp.setGroupingUsed(false);
        dfexp.setMaximumFractionDigits(4);

        DecimalFormat df = new DecimalFormat();
        df.setGroupingUsed(false);
        df.setMaximumFractionDigits(4);

        int cols = colNames.length;

        String[] str = new String[cols];

        Object[] obj = (Object[]) contents.get(row);
        int objlen = obj.length;
        for (int col = 0; col < cols; col++) {
            if (col >= objlen) {
                str[col] = "";
            } else if (obj[col] instanceof Double) {
                Double dbl = (Double) obj[col];
                if (dbl.doubleValue() >= .001 || dbl.doubleValue() == 0) {
                    str[col] = df.format(dbl.doubleValue());
                } else {
                    str[col] = dfexp.format(dbl.doubleValue());
                }
            } else {
                if (obj[col] == null) {
                    str[col] = "null";
                } else {
                    str[col] = obj[col].toString();
                }
            }
            if (str[col].length() > 9) {
                str[col].substring(0, 9);
            }
        }

        return str;

    }

    public String getTableTitle() {
        return title;
    }

    public String toString() {
        String out = "";
        String[][] str = (String[][]) getTableData();

        for (int i = 0; i < colNames.length; i++) {
            if (i > 0) {
                out += "\t";
            }
            String s = colNames[i] + "          ";
            out += s.substring(0, 10);
        }
        for (int i = 0; i < str.length; i++) {
            out += "\n";
            for (int j = 0; j < str[0].length; j++) {
                if (j > 0) {
                    out += "\t";
                }
                String s = str[i][j] + "          ";
                out += s.substring(0, 10);
            }
        }
        return out;
    }

    public List getContents() {
        return contents;
    }

    public int getNumberOfRows() {
        return contents.size();
    }

    public Object getValueAt(int row, int column) {
        return ((Object[]) contents.get(row))[column];
    }

    public void setValueAt(int row, int column, Object value) {
        ((Object[]) contents.get(row))[column] = value;
    }

    public int getRowCount() {
    	return contents.size();
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public int getColumnCount() {
        return colNames.length;
    }
}
