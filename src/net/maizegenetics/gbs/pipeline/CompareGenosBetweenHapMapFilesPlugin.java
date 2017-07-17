/*
 * CompareGenosBetweenHapMapFilesPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;

/**
 *
 * @author glaubitz
 */
public class CompareGenosBetweenHapMapFilesPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(CompareGenosBetweenHapMapFilesPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String hmp1FileStr, hmp2FileStr;
    private int startChr, endChr, chr, position;
    private HashMap<String, List<String>> taxaSynonyms = new HashMap<String, List<String>>();
    private HashMap<Integer, List<Integer>> taxaRedirect = new HashMap<Integer, List<Integer>>();
    private final String DELIMITER = "\t";

    public static enum SiteCompareType {

        SAME_STRAND, EITHER_STRAND, DIFF_STRAND, DIFFERENT
    };
    int nCompared = 0;
    int nSamePosNotComparable;
    // these are indices of the stat in the int[] array compareStats[]
    static final int NUM_TAXA_POSSIBLE_COMPARISONS = 0, NUM_TAXA_MISSING = 1, NUM_TAXA_COMPARED = 2, NUM_TAXA_DIFFERENT = 3, NUM_TAXA_HOMOZYGOUS_COMPARED = 4, NUM_TAXA_HOMOZYGOUS_DIFF = 5;
    //Taxa1:1	Taxa2:1	 NumberTested	Same	Different	ErrorRate	HomozygousTested	Same	Different	ErrorRate
    static final int COMPARE_STATS_LENGTH = 6; // this is the length of the int[] array compareStats[]
    static final int NUM_SITES_COMPARED = 0, NUM_SITES_DIFF = 1, NUM_SITES_HOMOZYGOUS_COMPARED = 2, NUM_SITES_HOMOZYGOUS_DIFF = 3;
    static final int COMPARE_TAXA_STATS_LENGTH = 4;
    static final int MINOR_ALLELE_FREQ1 = 0, MINOR_ALLELE_FREQ2 = 1, F_VALUE1 = 2, F_VALUE2 = 3; // these are the indices of the (double) summary stat in the double[] array summStats
    static final int summStatsLength = 4; //
    File outfile = null;
    DataOutputStream fw = null;
    private int myNumCalculations = 0;
    private List<Integer> myComparisons = new ArrayList<Integer>();
    private List<Double> myErrorRates = new ArrayList<Double>();
    private List<Integer> myHomComparisons = new ArrayList<Integer>();
    private List<Integer> myHomDiff = new ArrayList<Integer>();
    private List<Double> myHomError = new ArrayList<Double>();
    int[][][] myCompareStatsTaxa;

    public CompareGenosBetweenHapMapFilesPlugin() {
        super(null, false);
    }

    public CompareGenosBetweenHapMapFilesPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
                "\n\nThe options for CompareGenosBetweenHapMapFilesPlugin are:\n"
                + "    -hmp1  First hapmap format genotypic input file (use \"+\" as a wildcard character in place of the chromosome number)\n"
                + "    -hmp2  Second hapmap format genotypic input file to compare the first one to (use \"+\" as a wildcard character in place of the chromosome number)\n"
                + "    -sC    Start chromosome\n"
                + "    -eC    End chromosome\n"
                + "    -syn   Lookup table file of synonymous full taxon names in hmp1 and hmp2 (header line is ignored)\n"
                + "    -o     Output file for report (tab-delimited text format)\n\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        myLogger.addAppender(new ConsoleAppender(new SimpleLayout()));
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-hmp1", "--hapmap-file1", true);
            myArgsEngine.add("-hmp2", "--hapmap-file2", true);
            myArgsEngine.add("-sC", "--startChrom", true);
            myArgsEngine.add("-eC", "--endChrom", true);
            myArgsEngine.add("-syn", "--synonym-file", true);
            myArgsEngine.add("-o", "--output-file", true);
        }
        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-syn")) {
            String synFileStr = myArgsEngine.getString("-syn");
            File synFile = new File(synFileStr);
            if (!synFile.exists() || !synFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the file containing the synonym table for full taxon names (-syn option: " + synFileStr + ").");
            }
            if (!readTaxaSynonymsFromFile(synFile)) {
                throw new IllegalArgumentException("Problem reading the file containing the synonym table for full taxon names. Progam aborted.");
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a file containing the synonym table for full taxon names (option -syn).");
        }
        if (myArgsEngine.getBoolean("-sC")) {
            startChr = Integer.parseInt(myArgsEngine.getString("-sC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide a start chromosome (integer).\n");
        }
        if (myArgsEngine.getBoolean("-eC")) {
            endChr = Integer.parseInt(myArgsEngine.getString("-eC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide an end chromosome (integer).\n");
        }
        if (myArgsEngine.getBoolean("-hmp1")) {
            hmp1FileStr = myArgsEngine.getString("-hmp1");
            if ((hmp1FileStr.contains(File.separator) && hmp1FileStr.substring(hmp1FileStr.lastIndexOf(File.separator)).contains("+"))
                    || hmp1FileStr.contains("+")) {
                for (int chr = startChr; chr <= endChr; chr++) {
                    String infile = hmp1FileStr.replace("+", "" + chr);
                    File hmp1File = new File(infile);
                    if (!hmp1File.exists() || !hmp1File.isFile()) {
                        printUsage();
                        throw new IllegalArgumentException("Can't find the first hapmap format genotype input file (-hmp1 option: " + infile + ").");
                    }
                }
            } else {
                printUsage();
                throw new IllegalArgumentException("The name of the first hapmap input file should contain a \"+\" wildcard character in place of the chromosome number (-hmp1 option: " + hmp1FileStr + ").");
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the first hapmap format genotype input file (option -hmp1).");
        }
        if (myArgsEngine.getBoolean("-hmp2")) {
            hmp2FileStr = myArgsEngine.getString("-hmp2");
            if ((hmp2FileStr.contains(File.separator) && hmp2FileStr.substring(hmp2FileStr.lastIndexOf(File.separator)).contains("+"))
                    || hmp2FileStr.contains("+")) {
                for (int chr = startChr; chr <= endChr; chr++) {
                    String infile = hmp2FileStr.replace("+", "" + chr);
                    File hmp2File = new File(infile);
                    if (!hmp2File.exists() || !hmp2File.isFile()) {
                        printUsage();
                        throw new IllegalArgumentException("Can't find the second hapmap format genotype input file (-hmp2 option: " + infile + ").");
                    }
                }
            } else {
                printUsage();
                throw new IllegalArgumentException("The name of the second hapmap input file should contain a \"+\" wildcard character in place of the chromosome number (-hmp2 option: " + hmp2FileStr + ").");
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the second hapmap format genotype input file (option -hmp2).");
        }
        if (myArgsEngine.getBoolean("-o")) {
            String outFileStr = myArgsEngine.getString("-o");
            outfile = new File(outFileStr);
            try {
                fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));
                fw.writeBytes(
                        "Chr\tPosition\tAlleles1\tAlleles2\tCompareType\tMAF1\tMAF2\tf1\tf2\tSynonymousTaxaPairs\tMissing\tComparisons\tDiff\tErrorRate\tHomComparisons\tHomDiff\tHomErrorRate\n");
            } catch (Exception e) {
                throw new IllegalArgumentException("Unable to write to your output report file: " + e);
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output report file (inside an existing directory) (option -o).");
        }
    }

    @Override
    public DataSet performFunction(DataSet input) {
        for (chr = startChr; chr <= endChr; chr++) {
            String infile1 = hmp1FileStr.replace("+", "" + chr);
            String infile2 = hmp2FileStr.replace("+", "" + chr);
            myLogger.info("Comparing: " + infile1 + " to " + infile2);
            Alignment a1, a2;
            try {
                a1 = ImportUtils.readFromHapmap(infile1, null);
            } catch (Exception e) {
                myLogger.info("Could not read the first input hapmap file for chr" + chr + ":\n\t" + infile1 + "\n\tSkipping...");
                continue;
            }
            try {
                a2 = ImportUtils.readFromHapmap(infile2, null);
            } catch (Exception e) {
                myLogger.info("Could not read the second input hapmap file for chr" + chr + ":\n\t" + infile2 + "\n\tSkipping...");
                continue;
            }
            populateTaxaRedirect(a1, a2);
            myCompareStatsTaxa = new int[taxaRedirect.size()][][];
            findCommonPositionsAndCompare(a1, a2);
        }

        int[] comparisons = new int[myNumCalculations];
        double comparisonMean = 0.0;
        double[] errorRates = new double[myNumCalculations];
        double errorRateMean = 0.0;
        int[] homComparisons = new int[myNumCalculations];
        double homDiffSum = 0.0;
        double homComparisonSum = 0.0;
        double[] homErrors = new double[myNumCalculations];
        double homErrorSum = 0.0;
        for (int i = 0; i < myNumCalculations; i++) {
            comparisons[i] = myComparisons.get(i);
            comparisonMean = comparisonMean + comparisons[i];
            errorRates[i] = myErrorRates.get(i);
            errorRateMean = errorRateMean + errorRates[i];
            homComparisons[i] = myHomComparisons.get(i);
            homComparisonSum = homComparisonSum + homComparisons[i];
            homDiffSum = homDiffSum + myHomDiff.get(i);
            homErrors[i] = myHomError.get(i);
            homErrorSum = homErrorSum + homErrors[i];
        }

        double overAllHomoErrorRate = homDiffSum / homComparisonSum;

        comparisonMean = comparisonMean / myNumCalculations;
        double comparisonMedian = getMedian(comparisons);

        errorRateMean = errorRateMean / myNumCalculations;
        double errorRateMedian = getMedian(errorRates);

        double homComparisonMean = homComparisonSum / myNumCalculations;
        double homComparisonMedian = getMedian(homComparisons);

        double homErrorMean = homErrorSum / myNumCalculations;
        double homErrorMedian = getMedian(homErrors);

        myLogger.info("Comparison Mean\tComparison Median\tError Rate Mean\tError Rate Median\tHomozygous Comparison Mean\tHomozygous Comparison Median\tHomozygous Error Mean\tHomozygous Error Median");
        myLogger.info(comparisonMean + "\t" + comparisonMedian + "\t" + errorRateMean + "\t" + errorRateMedian + "\t" + homComparisonMean + "\t" + homComparisonMedian + "\t" + homErrorMean + "\t" + homErrorMedian);

        myLogger.info("Output Filename\tOver All Homo Error Rate\tCoverage");
        myLogger.info(outfile.getName() + "\t" + overAllHomoErrorRate + "\t" + (nSamePosNotComparable + nCompared));

        closeOutputFile();
        return null;
    }

    private static double getMedian(int[] values) {
        Arrays.sort(values);
        int middle = values.length / 2;
        if (values.length % 2 == 1) {
            return (values[middle - 1] + values[middle]) / 2.0;
        } else {
            return values[middle];
        }
    }

    private static double getMedian(double[] values) {
        Arrays.sort(values);
        int middle = values.length / 2;
        if (values.length % 2 == 1) {
            return (values[middle - 1] + values[middle]) / 2.0;
        } else {
            return values[middle];
        }
    }

    private boolean readTaxaSynonymsFromFile(File synFile) {
        taxaSynonyms.clear();
        String inputLine = "Nothing has been read from the taxon synonym input file yet";
        int nTaxa = 0, nCompareTaxa = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(synFile), 65536);
            inputLine = br.readLine();  // header line (ignore)
            while ((inputLine = br.readLine()) != null) {
                String[] cells = inputLine.split("\t");
                if (!(cells[0].equals("NA") || cells[1].equals("NA"))) {
                    List<String> synTaxaForTaxon = taxaSynonyms.get(cells[0]);
                    if (synTaxaForTaxon == null) {
                        taxaSynonyms.put(cells[0], synTaxaForTaxon = new ArrayList<String>());
                    }
                    synTaxaForTaxon.add(cells[1]);
                    ++nCompareTaxa;
                }
                ++nTaxa;
            }
        } catch (Exception e) {
            myLogger.error("Catch in reading taxon synonym input file e=" + e);
            e.printStackTrace();
            myLogger.info(inputLine);
            return false;
        }
        myLogger.info(nTaxa + " pairs of taxa full names read from the taxon synonym input file: " + nCompareTaxa + " of these will be compared (if found in the genotype files)");
        return true;
    }

    private void populateTaxaRedirect(Alignment a1, Alignment a2) {
        taxaRedirect.clear();
        myLogger.info("\n\nMaking list of comparable taxa found in the two hapmap files...\n");
        System.out.println("Taxon1\t\tTaxon2");
        System.out.println("------\t\t------");
        int nTaxaPairs = 0;
        for (int taxon1Index = 0; taxon1Index < a1.getSequenceCount(); taxon1Index++) {
            String taxon1 = a1.getFullTaxaName(taxon1Index);
            if (taxaSynonyms.containsKey(taxon1)) {
                for (String taxon2 : taxaSynonyms.get(taxon1)) {
                    for (int taxon2Index = 0; taxon2Index < a2.getSequenceCount(); taxon2Index++) {
                        if (taxon2.equals(a2.getFullTaxaName(taxon2Index))) {
                            List<Integer> synTaxaIndicesForTaxonIndex = taxaRedirect.get(taxon1Index);
                            if (synTaxaIndicesForTaxonIndex == null) {
                                taxaRedirect.put(taxon1Index, synTaxaIndicesForTaxonIndex = new ArrayList<Integer>());
                            }
                            synTaxaIndicesForTaxonIndex.add(taxon2Index);
                            System.out.println(taxon1 + "\t" + taxon2);
                            ++nTaxaPairs;
                            break;
                        }
                    }
                }
            }
        }
        myLogger.info("\nHapMap format genotype file1 contains " + a1.getSequenceCount() + " taxa in total\n");
        myLogger.info("\nHapMap format genotype file2 contains " + a2.getSequenceCount() + " taxa in total\n");
        myLogger.info("\n" + nTaxaPairs + " pairs of comparable taxa found in the two hapmap files\n\n");
    }

    private void findCommonPositionsAndCompare(Alignment a1, Alignment a2) {

        if (a1.getLoci().length != 1 || a2.getLoci().length != 1) {
            myLogger.error("ERROR: both hapmap genotype files should contain only a single chromosome");
            return;
        }
        if (!a1.getLoci()[0].getChromosomeName().equals(a2.getLoci()[0].getChromosomeName())) {
            myLogger.error("ERROR: the hapmap genotype files to compare do not contain the same chromosome");
            return;
        }
        if (Integer.parseInt(a1.getLoci()[0].getChromosomeName()) != chr || Integer.parseInt(a2.getLoci()[0].getChromosomeName()) != chr) {
            myLogger.error("ERROR: one or both of the hapmap genotype files to compare do not contain the expected chromosome "
                    + "(expected:" + chr + "  hmp1:" + a1.getLoci()[0].getChromosomeName() + "  hmp2:" + a2.getLoci()[0].getChromosomeName() + ")");
            return;
        }

        myLogger.info("\nHapMap format genotype file1 contains " + a1.getLocusSiteCount(a1.getLocus(0)) + " sites on chromosome " + a1.getLocusName(0) + "\n");
        myLogger.info("\nHapMap format genotype file2 contains " + a2.getLocusSiteCount(a2.getLocus(0)) + " sites on chromosome " + a2.getLocusName(0) + "\n\n");

        int nSites1 = a1.getSiteCount(), nSites2 = a2.getSiteCount();
        int s1 = 0, s2 = 0;
        nCompared = 0;
        nSamePosNotComparable = 0;
        boolean finished = false;
        while (!finished) {
            while (a1.getPositionInLocus(s1) != a2.getPositionInLocus(s2) && s1 < nSites1 && s2 < nSites2) {
                while (a1.getPositionInLocus(s1) < a2.getPositionInLocus(s2) && s1 < nSites1) {
                    s1++;
                }
                while (a2.getPositionInLocus(s2) < a1.getPositionInLocus(s1) && s2 < nSites2) {
                    s2++;
                }
            }
            if (s1 == nSites1 || s2 == nSites2) {
                finished = true;
            } else if (a1.getPositionInLocus(s1) == a2.getPositionInLocus(s2)) {
                position = a1.getPositionInLocus(s1);
                nCompared += getCompareTypeAndCompare(s1, a1, s2, a2);
                s2++;  // assumes that a1 is from GBS where duplicate SNPs have been removed (but a2 might contain some duplicates: e.g., hapmap2)
            }
        }
        myLogger.info(nCompared + " sites compared on chromosome " + chr
                + "\nAn additional " + nSamePosNotComparable + " sites on chromosome " + chr + " had the same position but incomparable alleles\n");

        outputTaxaReport(a1, a2);
    }

    private void outputTaxaReport(Alignment a1, Alignment a2) {

        System.out.println("Taxon1\tTaxon2\tNum_Sites_Compared\tNum_Sites_Diff\tNum_Sites_Homo_Compared\tNum_Sites_Homo_Diff");
        int taxon1Count = 0;
        for (Integer taxon1Index : taxaRedirect.keySet()) {
            List<Integer> synTaxaIndicesForTaxonIndex = taxaRedirect.get(taxon1Index);
            int taxon2Count = 0;
            for (Integer taxon2Index : synTaxaIndicesForTaxonIndex) {
                StringBuilder builder = new StringBuilder();
                builder.append(a1.getFullTaxaName(taxon1Index));
                builder.append(DELIMITER);
                builder.append(a2.getFullTaxaName(taxon2Index));
                builder.append(DELIMITER);
                builder.append(myCompareStatsTaxa[taxon1Count][taxon2Count][NUM_SITES_COMPARED]);
                builder.append(DELIMITER);
                builder.append(myCompareStatsTaxa[taxon1Count][taxon2Count][NUM_SITES_DIFF]);
                builder.append(DELIMITER);
                builder.append(myCompareStatsTaxa[taxon1Count][taxon2Count][NUM_SITES_HOMOZYGOUS_COMPARED]);
                builder.append(DELIMITER);
                builder.append(myCompareStatsTaxa[taxon1Count][taxon2Count][NUM_SITES_HOMOZYGOUS_DIFF]);
                System.out.println(builder.toString());
                taxon2Count++;
            }
            taxon1Count++;
        }

    }

    private int getCompareTypeAndCompare(int site1, Alignment a1, int site2, Alignment a2) {
        byte[] alleles1 = a1.getAlleles(site1);
        byte[] alleles2 = a2.getAlleles(site2);
        SiteCompareType compareType = getSiteCompareType(alleles1, alleles2);
        if (compareType == SiteCompareType.DIFFERENT) {
            nSamePosNotComparable++;
            return 0;
        }

        double[] summStats = new double[summStatsLength];
        summStats[MINOR_ALLELE_FREQ1] = a1.getMinorAlleleFrequency(site1);
        summStats[MINOR_ALLELE_FREQ2] = a2.getMinorAlleleFrequency(site2);
        summStats[F_VALUE1] = calculateF(a1, site1);
        summStats[F_VALUE2] = calculateF(a2, site2);
        String alleleString1 = a1.getBaseAsString(site1, alleles1[0]) + "/" + a1.getBaseAsString(site1, alleles1[1]);
        String alleleString2 = a2.getBaseAsString(site2, alleles2[0]) + "/" + a2.getBaseAsString(site2, alleles2[1]);

        int[][][] compareTaxaStatsSame = null;
        int[][][] compareTaxaStatsDiff = null;
        if ((compareType == SiteCompareType.SAME_STRAND) || (compareType == SiteCompareType.EITHER_STRAND)) {
            compareTaxaStatsSame = new int[taxaRedirect.size()][][];
        }

        if ((compareType == SiteCompareType.DIFF_STRAND) || (compareType == SiteCompareType.EITHER_STRAND)) {
            compareTaxaStatsDiff = new int[taxaRedirect.size()][][];
        }

        int[] compareStatsSame = null;
        int[] compareStatsDiff = null;
        int taxon1Count = 0;
        for (Integer taxon1Index : taxaRedirect.keySet()) {
            List<Integer> synTaxaIndicesForTaxonIndex = taxaRedirect.get(taxon1Index);
            if (myCompareStatsTaxa[taxon1Count] == null) {
                myCompareStatsTaxa[taxon1Count] = new int[synTaxaIndicesForTaxonIndex.size()][COMPARE_TAXA_STATS_LENGTH];
            }

            if ((compareType == SiteCompareType.SAME_STRAND) || (compareType == SiteCompareType.EITHER_STRAND)) {
                compareTaxaStatsSame[taxon1Count] = new int[synTaxaIndicesForTaxonIndex.size()][COMPARE_TAXA_STATS_LENGTH];
            }

            if ((compareType == SiteCompareType.DIFF_STRAND) || (compareType == SiteCompareType.EITHER_STRAND)) {
                compareTaxaStatsDiff[taxon1Count] = new int[synTaxaIndicesForTaxonIndex.size()][COMPARE_TAXA_STATS_LENGTH];
            }

            int taxon2Count = 0;
            for (Integer taxon2Index : synTaxaIndicesForTaxonIndex) {

                if ((compareType == SiteCompareType.SAME_STRAND) || (compareType == SiteCompareType.EITHER_STRAND)) {
                    int[] tempStats = compareGenotypes(taxon1Index, site1, a1, taxon2Index, site2, a2, true);
                    if (compareStatsSame == null) {
                        compareStatsSame = new int[COMPARE_STATS_LENGTH];
                    }
                    for (int i = 0; i < COMPARE_STATS_LENGTH; i++) {
                        compareStatsSame[i] += tempStats[i];
                    }
                    compareTaxaStatsSame[taxon1Count][taxon2Count][NUM_SITES_COMPARED] = 1;
                    compareTaxaStatsSame[taxon1Count][taxon2Count][NUM_SITES_DIFF] = tempStats[NUM_TAXA_DIFFERENT];
                    compareTaxaStatsSame[taxon1Count][taxon2Count][NUM_SITES_HOMOZYGOUS_COMPARED] = tempStats[NUM_TAXA_HOMOZYGOUS_COMPARED];
                    compareTaxaStatsSame[taxon1Count][taxon2Count][NUM_SITES_HOMOZYGOUS_DIFF] = tempStats[NUM_TAXA_HOMOZYGOUS_DIFF];
                }

                if ((compareType == SiteCompareType.DIFF_STRAND) || (compareType == SiteCompareType.EITHER_STRAND)) {
                    int[] tempStats = compareGenotypes(taxon1Index, site1, a1, taxon2Index, site2, a2, false);
                    if (compareStatsDiff == null) {
                        compareStatsDiff = new int[COMPARE_STATS_LENGTH];
                    }
                    for (int i = 0; i < COMPARE_STATS_LENGTH; i++) {
                        compareStatsDiff[i] += tempStats[i];
                    }
                    compareTaxaStatsDiff[taxon1Count][taxon2Count][NUM_SITES_COMPARED] = 1;
                    compareTaxaStatsDiff[taxon1Count][taxon2Count][NUM_SITES_DIFF] = tempStats[NUM_TAXA_DIFFERENT];
                    compareTaxaStatsDiff[taxon1Count][taxon2Count][NUM_SITES_HOMOZYGOUS_COMPARED] = tempStats[NUM_TAXA_HOMOZYGOUS_COMPARED];
                    compareTaxaStatsDiff[taxon1Count][taxon2Count][NUM_SITES_HOMOZYGOUS_DIFF] = tempStats[NUM_TAXA_HOMOZYGOUS_DIFF];
                }

                taxon2Count++;
            }
            taxon1Count++;
        }

        int[] compareResults;
        if (compareStatsSame == null) {
            compareResults = compareStatsDiff;
            addTaxaStats(myCompareStatsTaxa, compareTaxaStatsDiff);
        } else if (compareStatsDiff == null) {
            compareResults = compareStatsSame;
            addTaxaStats(myCompareStatsTaxa, compareTaxaStatsSame);
        } else {
            if (compareStatsSame[NUM_TAXA_DIFFERENT] <= compareStatsDiff[NUM_TAXA_DIFFERENT]) {
                compareResults = compareStatsSame;
                addTaxaStats(myCompareStatsTaxa, compareTaxaStatsSame);
            } else {
                compareResults = compareStatsDiff;
                addTaxaStats(myCompareStatsTaxa, compareTaxaStatsDiff);
            }
        }

        writeCompareStats(compareResults, alleleString1, alleleString2, compareType, summStats);

        return 1;
    }

    private void addTaxaStats(int[][][] totals, int[][][] additions) {
        for (int i = 0; i < totals.length; i++) {
            for (int j = 0; j < totals[0].length; j++) {
                for (int k = 0; k < totals[0][0].length; k++) {
                    totals[i][j][k] += additions[i][j][k];
                }
            }
        }
    }

    private double calculateF(Alignment a, int site) {
        byte majAllele = a.getMajorAllele(site);
        byte minAllele = a.getMinorAllele(site);
        int majGenoCnt = 0, minGenoCnt = 0, hetGenoCnt = 0;
        int nTaxa = a.getSequenceCount();
        // TERRY - Does this make sense?  What if it's het but not major/minor?
        for (int taxon = 0; taxon < nTaxa; taxon++) {
            byte[] bases = a.getBaseArray(taxon, site);
            if ((bases[0] == majAllele) && bases[1] == majAllele) {
                majGenoCnt++;
            } else if ((bases[0] == minAllele) && bases[1] == minAllele) {
                minGenoCnt++;
            } else if (AlignmentUtils.isEqual(bases, new byte[]{majAllele, minAllele})) {
                hetGenoCnt++;
            }
        }
        int nGenos = hetGenoCnt + majGenoCnt + minGenoCnt;
        if (nGenos > 0) {
            double propHets = (double) hetGenoCnt / nGenos;
            double maf = (double) (2 * majGenoCnt + hetGenoCnt) / (2 * nGenos);
            double expHets = 2.0 * maf * (1.0 - maf);
            return 1.0 - (propHets / expHets);
        } else {
            return Double.NaN;
        }
    }

    private int[] compareGenotypes(int taxon1Index, int site1, Alignment a1, int taxon2Index, int site2, Alignment a2, boolean sameStrand) {
        int[] compareStats = new int[COMPARE_STATS_LENGTH];
        byte base1 = a1.getBase(taxon1Index, site1);
        byte base2;
        if (sameStrand) {
            base2 = a2.getBase(taxon2Index, site2);
        } else {
            base2 = NucleotideAlignmentConstants.getNucleotideDiploidComplement(a2.getBase(taxon2Index, site2));
        }
        if (base1 != Alignment.UNKNOWN_DIPLOID_ALLELE && base2 != Alignment.UNKNOWN_DIPLOID_ALLELE) {
            if (!(AlignmentUtils.isHeterozygous(base1) || AlignmentUtils.isHeterozygous(base2))) {
                if (base1 != base2) {
                    ++compareStats[NUM_TAXA_DIFFERENT];
                    ++compareStats[NUM_TAXA_HOMOZYGOUS_DIFF];
                }
                ++compareStats[NUM_TAXA_HOMOZYGOUS_COMPARED];
            } else {
                if (!AlignmentUtils.isEqual(base1, base2)) {
                    ++compareStats[NUM_TAXA_DIFFERENT];
                }
            }
            ++compareStats[NUM_TAXA_COMPARED];
        } else {
            ++compareStats[NUM_TAXA_MISSING];
        }
        ++compareStats[NUM_TAXA_POSSIBLE_COMPARISONS];
        return compareStats;
    }

    private void writeCompareStats(int[] compareStats, String alleles1, String alleles2, SiteCompareType sct, double[] summStats) {
        double errRate = compareStats[NUM_TAXA_COMPARED] > 0 ? (double) compareStats[NUM_TAXA_DIFFERENT] / compareStats[NUM_TAXA_COMPARED] : Double.NaN;
        double errRateHom = compareStats[NUM_TAXA_HOMOZYGOUS_COMPARED] > 0 ? (double) compareStats[NUM_TAXA_HOMOZYGOUS_DIFF] / compareStats[NUM_TAXA_HOMOZYGOUS_COMPARED] : Double.NaN;

        try {
            fw.writeBytes(String.valueOf(chr));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(position));
            fw.writeBytes(DELIMITER);
            //fw.writeBytes((char) alleles1[0] + "/" + (char) alleles1[1] + "\t");
            //fw.writeBytes((char) alleles2[0] + "/" + (char) alleles2[1] + "\t");
            fw.writeBytes(alleles1);
            fw.writeBytes(DELIMITER);
            fw.writeBytes(alleles2);
            fw.writeBytes(DELIMITER);
            fw.writeBytes(sct.toString());
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(summStats[MINOR_ALLELE_FREQ1]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(summStats[MINOR_ALLELE_FREQ2]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(summStats[F_VALUE1]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(summStats[F_VALUE2]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(compareStats[NUM_TAXA_POSSIBLE_COMPARISONS]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(compareStats[NUM_TAXA_MISSING]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(compareStats[NUM_TAXA_COMPARED]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(compareStats[NUM_TAXA_DIFFERENT]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(errRate));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(compareStats[NUM_TAXA_HOMOZYGOUS_COMPARED]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(compareStats[NUM_TAXA_HOMOZYGOUS_DIFF]));
            fw.writeBytes(DELIMITER);
            fw.writeBytes(String.valueOf(errRateHom));
            fw.writeBytes("\n");

            myNumCalculations++;
            myComparisons.add(compareStats[NUM_TAXA_COMPARED]);
            myErrorRates.add(errRate);
            myHomDiff.add(compareStats[NUM_TAXA_HOMOZYGOUS_DIFF]);
            myHomComparisons.add(compareStats[NUM_TAXA_HOMOZYGOUS_COMPARED]);
            myHomError.add(errRateHom);

        } catch (Exception e) {
            throw new IllegalArgumentException("Unable to write to your output report file: " + e);
        }
    }

    private void closeOutputFile() {
        try {
            fw.close();
        } catch (Exception e) {
            throw new IllegalArgumentException("Unable to close your output report file: " + e);
        }
    }

    public static SiteCompareType getSiteCompareType(byte[] alleles1, byte[] alleles2) {

        //if (alleles1.length != 2 || alleles2.length != 2) {
        if (alleles1.length < 2 || alleles2.length < 2) {
            return SiteCompareType.DIFFERENT; // both must be biallelic to be compared
        }

        byte diploidValue1 = AlignmentUtils.getDiploidValue(alleles1[0], alleles1[1]);
        byte diploidValue2 = AlignmentUtils.getDiploidValue(alleles2[0], alleles2[1]);

        String iupac1 = NucleotideAlignmentConstants.getNucleotideIUPAC(diploidValue1);
        String iupac2 = NucleotideAlignmentConstants.getNucleotideIUPAC(diploidValue2);

        if ((iupac1 == null) || (iupac2 == null)) {
            return SiteCompareType.DIFFERENT;
        }

        //byte hetg1 = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(alleles1);
        //byte hetg2 = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(alleles2);

        char hetg1 = iupac1.charAt(0);
        char hetg2 = iupac2.charAt(0);

        if (hetg1 == '0' || hetg1 == '-' || hetg1 == '+' || hetg1 == Alignment.UNKNOWN_ALLELE_CHAR
                || hetg2 == '0' || hetg2 == '-' || hetg2 == '+' || hetg2 == Alignment.UNKNOWN_ALLELE_CHAR) {
            return SiteCompareType.DIFFERENT; // indel sites not compared
        }
        if (hetg1 == 'K' && hetg2 == 'K') {
            return SiteCompareType.SAME_STRAND;  // GT v GT
        }
        if (hetg1 == 'M' && hetg2 == 'M') {
            return SiteCompareType.SAME_STRAND;  // AC v AC
        }
        if (hetg1 == 'R' && hetg2 == 'R') {
            return SiteCompareType.SAME_STRAND;  // AG v AG
        }
        if (hetg1 == 'Y' && hetg2 == 'Y') {
            return SiteCompareType.SAME_STRAND;  // CT v CT
        }
        if (hetg1 == 'S' && hetg2 == 'S') {
            return SiteCompareType.EITHER_STRAND;  // CG v CG
        }
        if (hetg1 == 'W' && hetg2 == 'W') {
            return SiteCompareType.EITHER_STRAND;  // AT v AT
        }
        if (hetg1 == 'K' && hetg2 == 'M') {
            return SiteCompareType.DIFF_STRAND;  // GT v AC
        }
        if (hetg1 == 'M' && hetg2 == 'K') {
            return SiteCompareType.DIFF_STRAND;  // AC v GT
        }
        if (hetg1 == 'R' && hetg2 == 'Y') {
            return SiteCompareType.DIFF_STRAND;  // AG v CT
        }
        if (hetg1 == 'Y' && hetg2 == 'R') {
            return SiteCompareType.DIFF_STRAND;  // CT v AG
        }
        if (hetg1 == 'K' && hetg2 == 'R') {
            return SiteCompareType.DIFFERENT;  // GT v AG
        }
        if (hetg1 == 'K' && hetg2 == 'S') {
            return SiteCompareType.DIFFERENT;  // GT v CG
        }
        if (hetg1 == 'K' && hetg2 == 'W') {
            return SiteCompareType.DIFFERENT;  // GT v AT
        }
        if (hetg1 == 'K' && hetg2 == 'Y') {
            return SiteCompareType.DIFFERENT;  // GT v CT
        }
        if (hetg1 == 'M' && hetg2 == 'R') {
            return SiteCompareType.DIFFERENT;  // AC v AG
        }
        if (hetg1 == 'M' && hetg2 == 'S') {
            return SiteCompareType.DIFFERENT;  // AC v CG
        }
        if (hetg1 == 'M' && hetg2 == 'W') {
            return SiteCompareType.DIFFERENT;  // AC v AT
        }
        if (hetg1 == 'M' && hetg2 == 'Y') {
            return SiteCompareType.DIFFERENT;  // AC v CT
        }
        if (hetg1 == 'R' && hetg2 == 'K') {
            return SiteCompareType.DIFFERENT;  // AG v GT
        }
        if (hetg1 == 'R' && hetg2 == 'M') {
            return SiteCompareType.DIFFERENT;  // AG v AC
        }
        if (hetg1 == 'R' && hetg2 == 'S') {
            return SiteCompareType.DIFFERENT;  // AG v CG
        }
        if (hetg1 == 'R' && hetg2 == 'W') {
            return SiteCompareType.DIFFERENT;  // AG v AT
        }
        if (hetg1 == 'S' && hetg2 == 'K') {
            return SiteCompareType.DIFFERENT;  // CG v GT
        }
        if (hetg1 == 'S' && hetg2 == 'M') {
            return SiteCompareType.DIFFERENT;  // CG v AC
        }
        if (hetg1 == 'S' && hetg2 == 'R') {
            return SiteCompareType.DIFFERENT;  // CG v AG
        }
        if (hetg1 == 'S' && hetg2 == 'W') {
            return SiteCompareType.DIFFERENT;  // CG v AT
        }
        if (hetg1 == 'S' && hetg2 == 'Y') {
            return SiteCompareType.DIFFERENT;  // CG v CT
        }
        if (hetg1 == 'W' && hetg2 == 'K') {
            return SiteCompareType.DIFFERENT;  // AT v GT
        }
        if (hetg1 == 'W' && hetg2 == 'M') {
            return SiteCompareType.DIFFERENT;  // AT v AC
        }
        if (hetg1 == 'W' && hetg2 == 'R') {
            return SiteCompareType.DIFFERENT;  // AT v AG
        }
        if (hetg1 == 'W' && hetg2 == 'S') {
            return SiteCompareType.DIFFERENT;  // AT v CG
        }
        if (hetg1 == 'W' && hetg2 == 'Y') {
            return SiteCompareType.DIFFERENT;  // AT v CT
        }
        if (hetg1 == 'Y' && hetg2 == 'K') {
            return SiteCompareType.DIFFERENT;  // CT v GT
        }
        if (hetg1 == 'Y' && hetg2 == 'M') {
            return SiteCompareType.DIFFERENT;  // CT v AC
        }
        if (hetg1 == 'Y' && hetg2 == 'S') {
            return SiteCompareType.DIFFERENT;  // CT v CG
        }
        if (hetg1 == 'Y' && hetg2 == 'W') {
            return SiteCompareType.DIFFERENT;  // CT v AT
        }
        if (hetg1 == 'A' || hetg1 == 'C' || hetg1 == 'G' || hetg1 == 'T' || hetg2 == 'A' || hetg2 == 'C' || hetg2 == 'G' || hetg2 == 'T') {
            return SiteCompareType.DIFFERENT; // both must be variable to be compared (this is unlikely because then alleles.length==1)
        }

        return SiteCompareType.DIFFERENT; // default return
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
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
