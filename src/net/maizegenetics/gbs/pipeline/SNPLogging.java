/*
 * SNPLogging
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedWriter;
import java.io.File;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class SNPLogging {

    private final Logger myLogger = Logger.getLogger(SNPLogging.class);
    private static final String HEADER = "Chr\tPosition\tAlleles\tTagLocusStart\tStrand\tPlugin\tTest\tStatus\tValue\tCuttoff\n";
    private static final String DELIMITER = "\t";
    private final BufferedWriter myWriter;
    private final Class myCreator;

    public SNPLogging(String filename, Class creator) {

        if ((filename == null) || (filename.length() == 0)) {
            myWriter = null;
        } else {
            boolean exists = false;
            File file = new File(filename);
            if (file.exists()) {
                exists = true;
            }

            myWriter = Utils.getBufferedWriter(filename, true);
            if (!exists) {
                try {
                    myWriter.append(HEADER);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        myCreator = creator;

    }

    public void close() {
        try {
            myWriter.close();
        } catch (Exception e) {
            // do nothing;
        }
    }

    public void writeEntry(Alignment a, int site, String tagLocusStart, String strand, String test, String status, String value, String cuttoff) {
        String alleles = a.getMajorAlleleAsString(site) + "/" + a.getMinorAlleleAsString(site);
        writeEntry(a.getLocusName(site), a.getPositionInLocus(site), alleles, tagLocusStart, strand, test, status, value, cuttoff);
    }

    public void writeEntry(String chr, int position, String alleles, String tagLocusStart, String strand, String test, String status, String value, String cuttoff) {

        StringBuilder builder = new StringBuilder();
        builder.append(chr);
        builder.append(DELIMITER);
        builder.append(String.valueOf(position));
        builder.append(DELIMITER);
        builder.append(alleles);
        builder.append(DELIMITER);
        builder.append(tagLocusStart);
        builder.append(DELIMITER);
        builder.append(strand);
        builder.append(DELIMITER);
        builder.append(myCreator.getSimpleName());
        builder.append(DELIMITER);
        builder.append(test);
        builder.append(DELIMITER);
        builder.append(status);
        builder.append(DELIMITER);
        builder.append(value);
        builder.append(DELIMITER);
        builder.append(cuttoff);
        builder.append("\n");

        if (myWriter == null) {
            myLogger.info(builder.toString());
        } else {
            try {
                myWriter.append(builder.toString());
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }
}
