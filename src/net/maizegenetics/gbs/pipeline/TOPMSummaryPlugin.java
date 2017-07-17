/*
 * TOPMSummaryPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class TOPMSummaryPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(TOPMSummaryPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String myInputFilename = null;
    private TagsOnPhysicalMap myInputTOPM = null;
    private int myTagCount = 0;
    private int[] myChromosomes;
    private Map<Integer, Integer>[] myTagsPerSite;
    private Map<Integer, Set<Byte>>[] myVariantDefsPerPosition;
    private int myNumUndefinedStrandedTags = 0;
    private Set<Byte> myUndefinedStrandValues = new HashSet<Byte>();
    private String myOutputFilename = null;
    private int[] myNumTagsPerVariantsDefined;
    private TreeSet<Integer>[] myPositionsOnMaxVariantTags;

    public TOPMSummaryPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        myInputTOPM = new TagsOnPhysicalMap(myInputFilename, true);
        myTagCount = myInputTOPM.getTagCount();
        myLogger.info("performFunction: Number of Tags: " + myTagCount);

        myChromosomes = myInputTOPM.getChromosomes();
        Arrays.sort(myChromosomes);

        myNumTagsPerVariantsDefined = new int[myInputTOPM.getMaxNumVariants() + 1];

        myPositionsOnMaxVariantTags = new TreeSet[myChromosomes.length];
        for (int m = 0; m < myChromosomes.length; m++) {
            myPositionsOnMaxVariantTags[m] = new TreeSet<Integer>();
        }

        myTagsPerSite = new TreeMap[myChromosomes.length];
        for (int m = 0; m < myChromosomes.length; m++) {
            myTagsPerSite[m] = new TreeMap<Integer, Integer>();
        }
        myVariantDefsPerPosition = new TreeMap[myChromosomes.length];
        for (int m = 0; m < myChromosomes.length; m++) {
            myVariantDefsPerPosition[m] = new TreeMap<Integer, Set<Byte>>();
        }

        for (int i = 0; i < myTagCount; i++) {
            int startPos = myInputTOPM.getStartPosition(i);
            int endPos = myInputTOPM.getEndPosition(i);
            byte strand = myInputTOPM.getStrand(i);
            int chrom = myInputTOPM.getChromosome(i);
            int index = Arrays.binarySearch(myChromosomes, chrom);

            //String tag = BaseEncoder.getSequenceFromLong(myInputTOPM.getTag(i));
            if (strand == 1) {
                if (index < 0) {
                    myLogger.error("performFunction: tag: " + i + " chromosome: " + chrom + " not reported by getChromosomes()");
                    continue;
                }
                if (startPos > endPos) {
                    myLogger.error("performFunction: tag: " + i + " invalid state: strand: " + strand + "  start position: " + startPos + "  end position: " + endPos);
                    continue;
                }
                List<Integer> positionsOnTag = new ArrayList<Integer>();
                int numDefinedVariants = 0;
                for (int j = 0; j < myInputTOPM.getMaxNumVariants(); j++) {
                    int offset = myInputTOPM.getVariantPosOff(i, j);
                    byte def = myInputTOPM.getVariantDef(i, j);
                    if ((offset != Byte.MIN_VALUE) && (def != Byte.MIN_VALUE)) {
                        numDefinedVariants++;
                        int position = startPos + offset;
                        positionsOnTag.add(position);
                        Integer count = myTagsPerSite[index].get(position);
                        if (count == null) {
                            myTagsPerSite[index].put(position, 1);
                            Set<Byte> temp = new HashSet<Byte>();
                            temp.add(def);
                            myVariantDefsPerPosition[index].put(position, temp);
                        } else {
                            myTagsPerSite[index].put(position, count + 1);
                            Set<Byte> temp = myVariantDefsPerPosition[index].get(position);
                            temp.add(def);
                        }
                    }
                }
                if (numDefinedVariants == myInputTOPM.getMaxNumVariants()) {
                    myPositionsOnMaxVariantTags[index].addAll(positionsOnTag);
                }
                myNumTagsPerVariantsDefined[numDefinedVariants]++;
            } else if (strand == -1) {
                if (index < 0) {
                    myLogger.error("performFunction: tag: " + i + " chromosome: " + chrom + " not reported by getChromosomes()");
                    continue;
                }
                if (startPos < endPos) {
                    myLogger.error("performFunction: tag: " + i + " invalid state: strand: " + strand + "  start position: " + startPos + "  end position: " + endPos);
                    continue;
                }
                List<Integer> positionsOnTag = new ArrayList<Integer>();
                int numDefinedVariants = 0;
                for (int j = 0; j < myInputTOPM.getMaxNumVariants(); j++) {
                    int offset = myInputTOPM.getVariantPosOff(i, j);
                    byte def = myInputTOPM.getVariantDef(i, j);
                    if ((offset != Byte.MIN_VALUE) && (def != Byte.MIN_VALUE)) {
                        numDefinedVariants++;
                        int position = startPos + offset;
                        positionsOnTag.add(position);
                        Integer count = myTagsPerSite[index].get(position);
                        if (count == null) {
                            myTagsPerSite[index].put(position, 1);
                            Set<Byte> temp = new HashSet<Byte>();
                            temp.add(def);
                            myVariantDefsPerPosition[index].put(position, temp);
                        } else {
                            myTagsPerSite[index].put(position, count + 1);
                            Set<Byte> temp = myVariantDefsPerPosition[index].get(position);
                            temp.add(def);
                        }
                    }
                }
                if (numDefinedVariants == myInputTOPM.getMaxNumVariants()) {
                    myPositionsOnMaxVariantTags[index].addAll(positionsOnTag);
                }
                myNumTagsPerVariantsDefined[numDefinedVariants]++;
            } else {
                myNumUndefinedStrandedTags++;
                myUndefinedStrandValues.add(strand);
            }
        }

        for (int i = 0; i < myChromosomes.length; i++) {
            Iterator itr = myPositionsOnMaxVariantTags[i].iterator();
            StringBuilder builder = new StringBuilder();
            builder.append("performFunction: Chromosome: ");
            builder.append(myChromosomes[i]);
            builder.append(" Positions on Tags with Max Variants: ");
            boolean first = true;
            while (itr.hasNext()) {
                if (!first) {
                    builder.append(", ");
                } else {
                    first = false;
                }
                builder.append(itr.next());
            }
            myLogger.info(builder.toString());
        }

        myLogger.info("performFunction: Number of Tags with Undefined Strands: " + myNumUndefinedStrandedTags);
        Iterator itr = myUndefinedStrandValues.iterator();
        while (itr.hasNext()) {
            myLogger.info("performFunction: Undefined Strand Value: " + itr.next());
        }

        int totalSNPs = 0;
        for (int i = 0; i < myChromosomes.length; i++) {
            totalSNPs += myTagsPerSite[i].size();
            myLogger.info("performFunction: Chromosome: " + myChromosomes[i] + " Number of SNPs: " + myTagsPerSite[i].size());
        }
        myLogger.info("performFunction: Total SNPs: " + totalSNPs);

        for (int i = 0; i <= myInputTOPM.getMaxNumVariants(); i++) {
            myLogger.info("performFunction: Number of Tags: " + myNumTagsPerVariantsDefined[i] + " Has: " + i + " Variants Defined");
        }

        printSummary();
        return null;
    }

    private void printSummary() {
        BufferedWriter writer = null;

        try {
            writer = Utils.getBufferedWriter(myOutputFilename);
            writer.append("Chromosome\tPosition\tNum Tags\tVariant Defs\n");
            for (int c = 0; c < myChromosomes.length; c++) {
                Iterator itr = myTagsPerSite[c].entrySet().iterator();
                while (itr.hasNext()) {
                    Map.Entry entry = (Map.Entry) itr.next();
                    writer.append(myChromosomes[c] + "\t" + entry.getKey() + "\t" + entry.getValue() + "\t");
                    Set<Byte> defSet = myVariantDefsPerPosition[c].get(entry.getKey());
                    Iterator itr2 = defSet.iterator();
                    boolean notFirst = false;
                    while (itr2.hasNext()) {
                        if (notFirst) {
                            writer.append(",");
                        } else {
                            notFirst = true;
                        }
                        writer.append(NucleotideAlignmentConstants.getHaplotypeNucleotide(((Byte) itr2.next()).byteValue()));
                    }
                    writer.append("\n");
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                writer.close();
            } catch (Exception ex) {
                // do nothing
            }
        }
    }

    private void printUsage() {
        myLogger.info(
                "\nThe options for the TOPMSummaryPlugin:\n"
                + "-input Input TOPM\n"
                + "-output Output Filename\n");
    }

    @Override
    public void setParameters(String[] args) {

        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-input", "-input", true);
            myArgsEngine.add("-output", "-output", true);
        }
        myArgsEngine.parse(args);

        myInputFilename = myArgsEngine.getString("-input");
        if ((myInputFilename == null) || (myInputFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: Must define input file");
        }
        File inputFile = new File(myInputFilename);
        if (!inputFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: The input file doesn't exist: " + myInputFilename);
        }

        myOutputFilename = myArgsEngine.getString("-output");
        if ((myOutputFilename == null) || (myOutputFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: Must define output file");
        }
        File outputFile = new File(myOutputFilename);
        if (outputFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: The output file already exists: " + myOutputFilename);
        }

    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
