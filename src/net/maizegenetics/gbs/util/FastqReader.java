/*
 * FastqReader
 */
package net.maizegenetics.gbs.util;

import java.io.BufferedReader;
import net.maizegenetics.util.Utils;

/**
 *
 * @author terry
 */
public class FastqReader {
    
    private BufferedReader myReader;
    private String myFilename;
    
    public FastqReader(String filename) {
        myFilename = filename;
        myReader = Utils.getBufferedReader(myFilename);
    }

    /**
     * Returns single entry in fastq file consisting of 4 Strings.
     *
     * @return Array of four Strings (Null if end of file).
     */
    public String[] getNext() {

        // Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
        // Line 2 is the raw sequence letters.
        // Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
        // Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

        try {
            String[] result = new String[4];
            result[0] = myReader.readLine();
            if (result[0] == null) {
                return null;
            }
            result[1] = myReader.readLine();
            result[2] = myReader.readLine();
            result[3] = myReader.readLine();
            return result;
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        
    }
    
    public void close() {
        try {
            myReader.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
