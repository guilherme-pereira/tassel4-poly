/*
 * SeparatePlugin.java
 *
 * Created on July 31, 2010
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import java.awt.Frame;

import javax.swing.*;

import java.net.URL;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class SeparatePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(SeparatePlugin.class);
    private String[] myChromosomesToSeparate = null;

    /**
     * Creates a new instance of SeparatePlugin
     */
    public SeparatePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> inputs = input.getDataSet();
            List<DataSet> result = new ArrayList<DataSet>();

            for (Datum current : inputs) {
                Object currentValue = current.getData();
                if (currentValue instanceof MarkerPhenotype) {

                    MarkerPhenotype mp = (MarkerPhenotype) currentValue;
                    Phenotype pheno = mp.getPhenotype();
                    String phenoName = current.getName() + "_pheno";
                    Datum phenoDatum = new Datum(phenoName, pheno, null);

                    Alignment align = mp.getAlignment();
                    String alignName = current.getName() + "_align";
                    Datum alignDatum = new Datum(alignName, align, null);

                    DataSet tds = new DataSet(new Datum[]{phenoDatum, alignDatum}, this);
                    result.add(tds);

                } else if (currentValue instanceof Alignment) {

                    List<Datum> alignments = separateAlignmentIntoLoci((Alignment) currentValue, current.getName(), myChromosomesToSeparate);
                    if (alignments.size() > 0) {
                        DataSet tds = new DataSet(alignments, this);
                        result.add(tds);
                    }

                }
            }

            if (result.isEmpty()) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Nothing to Separate");
                } else {
                    myLogger.warn("performFunction: Nothing to Separate.");
                }
                return null;
            } else {
                DataSet resultDataSet = DataSet.getDataSet(result, this);
                fireDataSetReturned(new PluginEvent(resultDataSet, SeparatePlugin.class));
                return resultDataSet;
            }
        } finally {
            fireProgress(100);
        }

    }

    public static List<Datum> separateAlignmentIntoLoci(Alignment alignment, String dataSetName) {
        return separateAlignmentIntoLoci(alignment, dataSetName, null);
    }

    public static List<Datum> separateAlignmentIntoLoci(Alignment alignment, String dataSetName, String[] chromosomesToSeparate) {

        List<Datum> result = new ArrayList<Datum>();
        Alignment[] alignments = alignment.getAlignments();
        for (int i = 0; i < alignments.length; i++) {
            int[] offsets = alignments[i].getLociOffsets();
            if (offsets.length > 1) {
                Locus[] loci = alignments[i].getLoci();
                for (int j = 0; j < offsets.length; j++) {
                    if (alignmentInList(loci[j], chromosomesToSeparate)) {
                        String name;
                        if (dataSetName == null) {
                            name = "Alignment_chrom" + loci[j];
                        } else {
                            name = dataSetName + "_chrom" + loci[j];
                        }
                        int endSite;
                        try {
                            endSite = offsets[j + 1] - 1;
                        } catch (Exception e) {
                            endSite = alignments[i].getSiteCount() - 1;
                        }
                        Datum td = new Datum(name, FilterAlignment.getInstance(alignments[i], offsets[j], endSite), null);
                        result.add(td);
                    }
                }
            } else {
                if ((alignments.length > 1) && (alignmentInList(alignments[i].getLoci()[0], chromosomesToSeparate))) {
                    String name;
                    if (dataSetName == null) {
                        name = "Alignment_chrom" + alignments[i].getLocus(0);
                    } else {
                        name = dataSetName + "_chrom" + alignments[i].getLocus(0);
                    }
                    Datum td = new Datum(name, alignments[i], null);
                    result.add(td);
                }
            }
        }

        return result;

    }

    private static boolean alignmentInList(Locus locus, String[] chromosomesToSeparate) {

        if (chromosomesToSeparate == null) {
            return true;
        }

        String currentChr = locus.getChromosomeName();
        for (int i = 0; i < chromosomesToSeparate.length; i++) {
            if (currentChr.equalsIgnoreCase(chromosomesToSeparate[i])) {
                return true;
            }
        }

        return false;
    }

    public void setChromosomesToSeparate(String[] chrs) {
        myChromosomesToSeparate = new String[chrs.length];
        for (int i = 0; i < chrs.length; i++) {
            myChromosomesToSeparate[i] = chrs[i].trim();
        }
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = SeparatePlugin.class.getResource("images/Separate.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Separate";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Separate Data (i.e. into Chromosomes)";
    }
}
