/*
 * ConvertAlignmentCoordinatesPlugin
 */
package net.maizegenetics.baseplugins;

import java.awt.Frame;

import java.io.BufferedReader;

import java.util.HashMap;
import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import java.util.List;
import java.util.regex.Pattern;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MutableAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class ConvertAlignmentCoordinatesPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ConvertAlignmentCoordinatesPlugin.class);
    private String myMapFilename = null;
    private final HashMap<String, Locus> myLociMap = new HashMap<String, Locus>();
    private final HashMap<String, Locus> myAlignmentLociMap = new HashMap<String, Locus>();

    public ConvertAlignmentCoordinatesPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(Alignment.class);

            if (alignInList.size() != 1) {
                String gpMessage = "Invalid selection.  Please select one genotype alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), gpMessage);
                } else {
                    myLogger.error(gpMessage);
                }
                return null;
            }
            Datum current = alignInList.get(0);
            Datum result = processDatum(current);
            DataSet output = new DataSet(result, this);
            fireDataSetReturned(new PluginEvent(output, FilterAlignmentPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }

    }

    private Datum processDatum(Datum input) {

        if ((myMapFilename == null) || (myMapFilename.length() == 0)) {
            throw new IllegalStateException("ConvertAlignmentCoordinatesPlugin: processDatum: map filename not set.");
        }

        MutableAlignment alignment = null;

        if (input.getData() instanceof MutableAlignment) {
            alignment = (MutableAlignment) input.getData();
        } else {
            alignment = MutableNucleotideAlignment.getInstance((Alignment) input.getData());
        }

        Locus[] loci = alignment.getLoci();
        myLociMap.clear();
        for (int i = 0; i < loci.length; i++) {
            myLociMap.put(loci[i].getName(), loci[i]);
            myAlignmentLociMap.put(loci[i].getName(), loci[i]);
        }

        String[] snpIDs = alignment.getSNPIDs();

        BufferedReader br = null;
        int count = 1;
        try {
            br = Utils.getBufferedReader(myMapFilename);
            Pattern sep = Pattern.compile("\\s+");

            // Ignore head, then get first line.
            String inputline = br.readLine();
            inputline = br.readLine();
            int numChanges = 0;
            while (inputline != null) {

                count++;
                inputline = inputline.trim();
                // ID	v1_chr	v1_pos	v2_chr	v2_pos	strand
                String[] parsedline = sep.split(inputline);
                inputline = br.readLine();

                String locus1 = getLocusName(parsedline[1]);

                if (myAlignmentLociMap.get(locus1) != null) {

                    String snpID = new String(parsedline[0]);
                    int pos1 = Integer.valueOf(parsedline[2]);
                    String locus2 = getLocusName(parsedline[3]);
                    int pos2 = Integer.valueOf(parsedline[4]);

                    if ((!locus1.equals(locus2)) || (pos1 != pos2)) {

                        int site = getSiteOfSNPID(snpID, snpIDs);
                        if (site < 0) {
                            //myLogger.warn("map file line: " + count + "  SNP ID: " + snpID + "  position: " + pos1 + "  locus: " + locus1 + " not found in alignment.");
                            continue;
                        }

                        if ((pos1 != alignment.getPositionInLocus(site)) || (!locus1.equals(alignment.getLocus(site).getName()))) {
                            myLogger.warn("map file line: " + count + "  SNP ID: " + snpID + "  position: " + pos1 + "  locus: " + locus1 + " position and locus do not match alignment.");
                            myLogger.warn("Alignment SNP ID: " + alignment.getSNPID(site) + "  position: " + alignment.getPositionInLocus(site) + "  locus: " + alignment.getLocusName(site));
                            continue;
                        }

                        numChanges++;

                        alignment.setPositionOfSite(site, pos2);

                        if (!locus1.equals(locus2)) {
                            alignment.setLocusOfSite(site, getLocusObj(locus2));
                        }

                    }

                }

            }

            myLogger.info("Number Changes: " + numChanges);

            alignment.clean();

            return new Datum(input.getName() + "_NewCoordinates", alignment, null);

        } catch (Exception e) {
            myLogger.error("processDatum: problem converting alignment: line: " + count + "  message: " + e.getMessage());
        } finally {
            try {
                br.close();
            } catch (Exception e) {
                // do nothing
            }
        }

        return null;
    }

    private Locus getLocusObj(String locus) {
        Locus locusObj = myLociMap.get(locus);
        if (locusObj == null) {
            locusObj = new Locus(locus, locus, -1, -1, null, null);
            myLociMap.put(locus, locusObj);
        }
        return locusObj;
    }

    private String getLocusName(String input) {
        if (input.startsWith("chr")) {
            return input.substring(3);
        } else {
            return new String(input);
        }
    }

    private int getSiteOfSNPID(String id, String[] snpIDs) {
        for (int i = 0, n = snpIDs.length; i < n; i++) {
            if (id.equals(snpIDs[i])) {
                return i;
            }
        }
        return -1;
    }

    public void setMapFilename(String filename) {
        myMapFilename = filename;
    }

    public String getMapFilename() {
        return myMapFilename;
    }

    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
