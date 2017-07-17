package net.maizegenetics.stats.GLM;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Apr 29, 2004
 * Time: 9:35:35 AM
 *
 * Provides methods to get data required by a ReportWriter
 *
 */
public interface Reporter {
    public String[] getColumnNames();
    public Object[][] getResults();
}