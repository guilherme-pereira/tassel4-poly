/*
 * MergeAlignmentsSameSitesPlugin
 */
package net.maizegenetics.baseplugins;

import java.awt.Frame;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.ImageIcon;

import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class MergeAlignmentsSameSitesPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(MergeAlignmentsSameSitesPlugin.class);
    private static final char DELIMITER = '\t';
    private static final Pattern DELIMITER_PATTERN = Pattern.compile(String.valueOf(DELIMITER));
    private List<String> myInputFiles;
    private String myOutputFile;

    public MergeAlignmentsSameSitesPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public DataSet performFunction(DataSet input) {

        if ((myInputFiles == null) || (myInputFiles.size() < 2)) {
            myLogger.warn("performFunction: Must set at least two hapmap alignment files.");
            return null;
        }

        if ((myOutputFile == null) || (myOutputFile.length() == 0)) {
            myLogger.warn("performFunction: Must set an output file.");
            return null;
        }

        FileWriter fw = null;
        BufferedWriter bw = null;
        int numInputs = myInputFiles.size();
        BufferedReader[] readers = new BufferedReader[numInputs];
        try {

            fw = new FileWriter(Utils.addSuffixIfNeeded(myOutputFile, ".hmp.txt"));
            bw = new BufferedWriter(fw, 1000000);

            writeHeader(bw);
            for (int i = 0; i < numInputs; i++) {
                readers[i] = Utils.getBufferedReader(myInputFiles.get(i));
                String header = readers[i].readLine();
                int index = Utils.findNthOccurrenceInString(header, DELIMITER, ImportUtils.NUM_HAPMAP_NON_TAXA_HEADERS);
                bw.write(header.substring(index));
            }
            bw.write("\n");

            String[] current = new String[numInputs];
            int line = 1;
            while ((current[0] = readers[0].readLine()) != null) {
                line++;
                int firstIndex = Utils.findNthOccurrenceInString(current[0], DELIMITER, ImportUtils.NUM_HAPMAP_NON_TAXA_HEADERS);
                String[] firstHeader = DELIMITER_PATTERN.split(current[0].substring(0, firstIndex));
                bw.write(current[0]);
                for (int i = 1; i < numInputs; i++) {
                    current[i] = readers[i].readLine();
                    int currentIndex = Utils.findNthOccurrenceInString(current[i], DELIMITER, ImportUtils.NUM_HAPMAP_NON_TAXA_HEADERS);
                    String[] currentHeader = DELIMITER_PATTERN.split(current[i].substring(0, currentIndex));
                    if (!firstHeader[0].equals(currentHeader[0])) {
                        throw new IllegalStateException("MergeAlignmentsSameSitesPlugin: performFunction: Site Name does not match line: " + line + "  first: " + firstHeader[0] + "  current: " + currentHeader[0]);
                    }
                    if (!firstHeader[2].equals(currentHeader[2])) {
                        throw new IllegalStateException("MergeAlignmentsSameSitesPlugin: performFunction: Chromosome does not match line: " + line + "  first: " + firstHeader[2] + "  current: " + currentHeader[2]);
                    }
                    if (!firstHeader[3].equals(currentHeader[3])) {
                        throw new IllegalStateException("MergeAlignmentsSameSitesPlugin: performFunction: Physical Position does not match line: " + line + "  first: " + firstHeader[3] + "  current: " + currentHeader[3]);
                    }
                    bw.write(current[i].substring(currentIndex));
                }
                bw.write("\n");
            }


        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }

            for (int i = 0; i < readers.length; i++) {
                try {
                    readers[i].close();
                } catch (Exception e) {
                    // do nothing
                }
            }
            fireProgress(100);
        }

        return null;

    }

    private void writeHeader(BufferedWriter bw) {
        try {
            bw.write("rs#");
            bw.write(DELIMITER);
            bw.write("alleles");
            bw.write(DELIMITER);
            bw.write("chrom");
            bw.write(DELIMITER);
            bw.write("pos");
            bw.write(DELIMITER);
            bw.write("strand");
            bw.write(DELIMITER);
            bw.write("assembly#");
            bw.write(DELIMITER);
            bw.write("center");
            bw.write(DELIMITER);
            bw.write("protLSID");
            bw.write(DELIMITER);
            bw.write("assayLSID");
            bw.write(DELIMITER);
            bw.write("panelLSID");
            bw.write(DELIMITER);
            bw.write("QCcode");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void setInputFiles(List<String> files) {
        myInputFiles = new ArrayList<String>(files);
    }

    public void setOutputFile(String filename) {
        myOutputFile = filename;
    }

    public ImageIcon getIcon() {
        return null;
    }

    public String getButtonName() {
        return "Merge";
    }

    public String getToolTipText() {
        return "Merge Alignments Same Sites";
    }
}