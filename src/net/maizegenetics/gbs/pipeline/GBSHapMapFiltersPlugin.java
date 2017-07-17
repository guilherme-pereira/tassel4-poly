package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.File;

import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.ImageIcon;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.VCFUtil;

import org.apache.log4j.Logger;

/**
 * Basic filters needed for removing bad sites and taxa from GBS pipelines
 *
 * @author edbuckler
 */
public class GBSHapMapFiltersPlugin extends AbstractPlugin {

    private int startChromosome = 1, endChromosome = 10;
    private ArgsEngine myArgsEngine = null;
    private static final Logger myLogger = Logger.getLogger(GBSHapMapFiltersPlugin.class);
    private String snpLogFileName;
    private SNPLogging snpLogging = null;
    private String suppliedInputFileName, suppliedOutputFileName, infile, outfile;
    private double minF = -2.0, minMAF = 0, maxMAF = 1, minPresence = 0;
    private boolean usePedigree = false;
    private HashMap<String, Double> taxaFs = null;
    private boolean hLD = false;
    private double minR2 = 0.01;
    private double minBonP = 0.01;
    private String[] lowCoverageTaxa = null;
    private static enum INPUT_FORMAT {hapmap, vcf}; //input file format, acceptable values are "hapmap" "vcf" 
    private INPUT_FORMAT inputFormat = INPUT_FORMAT.hapmap;
    private int myMaxNumAlleles=VCFUtil.VCF_DEFAULT_MAX_NUM_ALLELES;

    public GBSHapMapFiltersPlugin() {
        super(null, false);
    }

    public GBSHapMapFiltersPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        for (int chr = startChromosome; chr <= endChromosome; chr++) {
            infile = suppliedInputFileName.replace("+", "" + chr);
            outfile = suppliedOutputFileName.replace("+", "" + chr);
            myLogger.info("Reading: " + infile);
            Alignment a;
            
            if (inputFormat == INPUT_FORMAT.hapmap)
            {
                try {
                    a = ImportUtils.readFromHapmap(infile, this);
                } catch (Exception e) {
                    myLogger.info("Could not read input hapmap file for chr" + chr + ":\n\t" + infile + "\n\tSkipping...");
                    continue;
                }
            }
            else if (inputFormat == INPUT_FORMAT.vcf)
            {
                try {
                    a = ImportUtils.readFromVCF(infile, this, myMaxNumAlleles);
                } catch (Exception e) {
                    myLogger.info("Could not read input vcf file for chr" + chr + ":\n\t" + infile + "\n\tSkipping...");
                    continue;
                }
            }
            else
            {
                 throw new IllegalArgumentException("File format " + inputFormat + " is not recognized!");
            }
            
            myLogger.info("Original Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            if (a.getSiteCount() == 0) {
                continue;
            }
            double realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, false, false);
            double randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, true, false);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(a, false);
            myLogger.info("Ratio of RandomToReal:" + randomDist / realDist);
            if (myArgsEngine.getBoolean("-mnTCov")) {
                double tCov = Double.parseDouble(myArgsEngine.getString("-mnTCov"));
                if (lowCoverageTaxa == null) {
                    lowCoverageTaxa = getLowCoverageLines(a, tCov);
                }  // Note: lowCoverageTaxa is based upon the startChromosome only
                IdGroup keepTaxa = AlignmentFilterByGBSUtils.getFilteredIdGroupByName(a.getIdGroup(), lowCoverageTaxa, false);
                a = FilterAlignment.getInstance(a, keepTaxa);
                myLogger.info("TaxaFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
                if (a.getSiteCount() == 0) {
                    continue;
                }
            }
            realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, false, false);
            randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, true, false);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(a, false);
            myLogger.info("Ratio of RandomToReal:" + randomDist / realDist);
            int minCount = (int) Math.round(a.getSequenceCount() * minPresence);
            if (usePedigree) {
                // filter the sites for minCount, minMAF and maxMAF (but not minF) based on all of the taxa
                int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, false, -2.0, minCount, minMAF, maxMAF, snpLogging, "Filter the sites for minCount, minMAF and maxMAF (but not minF) based on all of the taxa");
                a = FilterAlignment.getInstance(a, goodLowHetSites);

                // filter the sites for minF only based only on the taxa with expectedF >= minF
                String[] highExpectedFTaxa = getHighExpectedFTaxa(a);
                IdGroup highExpectedFTaxaIDGroup = AlignmentFilterByGBSUtils.getFilteredIdGroupByName(a.getIdGroup(), highExpectedFTaxa, true);
                Alignment inbredGenos = FilterAlignment.getInstance(a, highExpectedFTaxaIDGroup);
                int[] goodLowFSites = AlignmentFilterByGBSUtils.getLowHetSNPs(inbredGenos, false, minF, 0, -0.1, 2.0, snpLogging, "Filter the sites for minF only based only on the taxa with expectedF >= minF");
                inbredGenos = null;
                System.gc();
                a = FilterAlignment.getInstance(a, goodLowFSites);
            } else {
                int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, false, minF, minCount, minMAF, maxMAF, snpLogging, "Filter the sites");
                a = FilterAlignment.getInstance(a, goodLowHetSites);
            }
            myLogger.info("SiteFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            if (a.getSiteCount() == 0) {
                continue;
            }
            realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, false, true);
            randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, true, false);
            System.out.printf("%d %d %g %g %g %n", a.getSiteCount(), minCount, minF, minMAF, randomDist / realDist);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(a, false);
            
            if (inputFormat == INPUT_FORMAT.hapmap)
            {
                ExportUtils.writeToHapmap(a, false, outfile, '\t', null);
            }
            else
            {
                ExportUtils.writeToVCF(a, outfile, '\t');
            }
            
            myLogger.info("File written after basic filtering:" + outfile);

            if (hLD) {
                a = ImportUtils.readFromHapmap(outfile, null);
                int[] gs = AlignmentFilterByGBSUtils.getGoodSitesByLD(a, minR2, minBonP, 128, 100, 20, false);
                a = FilterAlignment.getInstance(a, gs);
                myLogger.info("LDFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
                if (a.getSiteCount() == 0) {
                    continue;
                }
                ExportUtils.writeToHapmap(a, false, outfile, '\t', null);
                myLogger.info("File written after basic & LD filtering:" + outfile);
            }
        }
        snpLogging.close();
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\n\n\nThe GBSHapMapFiltersPlugin accepts the following options:\n"
                + "-hmp     Input HapMap file; use a plus sign (+) as a wild card character to \n"
                + "           specify multiple chromosome numbers.\n"
                + "-vcf     Input VCF file. Use a plus sign (+) as a wild card character to specify \n"
                + "         multiple chromosome numbers. Options -hmp and -vcf are mutual exclusive.\n"
                + "-o       Output HapMap file\n"
                + "-mnTCov  Minimum taxa coverage (default: no filter)\n"
                + "-mnSCov  Minimum presence (default: no filter)\n"
                + "-mnF     Minimum F (inbreeding coefficient) (default -2.0 = no filter)\n"
                + "-p       Pedigree file containing full sample names (or expected names after merging) & expected inbreeding\n"
                + "         coefficient (F) for each.  Only taxa with expected F >= mnF used to calculate F = 1-Ho/He.\n"
                + "         (default: use ALL taxa to calculate F)\n"
                + "-mnMAF   Minimum minor allele frequency (default: 0.0 = no filter)\n"
                + "-mxMAF   Maximum minor allele frequency (default: 1.0 = no filter)\n"
                + "-hLD     Filter for high LD\n"
                + "-mnR2    Minimum R-square value for the LD filter (default: " + minR2 + ")\n"
                + "-mnBonP  Minimum Bonferroni-corrected p-value for the LD filter (default: " + minBonP + ")\n"
                + "-sC      Start chromosome (default: 1).\n"
                + "-eC      End chromosome (default: 10).\n"
                + "-maxAlleleVCF   Maximum number of alleles allowed in vcf file.\n"
                + "-snpLog  SNPs Removed Log file name\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-hmp", "-hmpFile", true);
            myArgsEngine.add("-vcf", "-vcfFile", true);
            myArgsEngine.add("-o", "--outFile", true);
            myArgsEngine.add("-mnTCov", "--minTaxaCov", true);
            myArgsEngine.add("-mnSCov", "--minSiteCov", true);
            myArgsEngine.add("-mnF", "--minFInbreeding", true);
            myArgsEngine.add("-p", "--pedigree-file", true);
            myArgsEngine.add("-mnMAF", "--minMinorAlleleFreq", true);
            myArgsEngine.add("-mxMAF", "--maxinorAlleleFreq", true);
            myArgsEngine.add("-hLD", "--highLD", false);
            myArgsEngine.add("-mnR2", "--minRSquare", true);
            myArgsEngine.add("-mnBonP", "--minBonferronPForLD", true);
            myArgsEngine.add("-sC", "--startChrom", true);
            myArgsEngine.add("-eC", "--endChrom", true);
            myArgsEngine.add("-maxAlleleVCF", "--maxAlleleVCF", true);
            myArgsEngine.add("-snpLog", "", true);
        }

        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-sC")) {
            startChromosome = Integer.parseInt(myArgsEngine.getString("-sC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide a start chromosome.\n");
        }
        if (myArgsEngine.getBoolean("-eC")) {
            endChromosome = Integer.parseInt(myArgsEngine.getString("-eC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide an end chromosome.\n");
        }
        if (myArgsEngine.getBoolean("-mnSCov")) {
            minPresence = Double.parseDouble(myArgsEngine.getString("-mnSCov"));
        }
        if (myArgsEngine.getBoolean("-mnF")) {
            minF = Double.parseDouble(myArgsEngine.getString("-mnF"));
        }
        if (myArgsEngine.getBoolean("-p")) {
            String pedigreeFileStr = myArgsEngine.getString("-p");
            File pedigreeFile = new File(pedigreeFileStr);
            if (!pedigreeFile.exists() || !pedigreeFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the pedigree input file (-p option: " + pedigreeFileStr + ").");
            }
            taxaFs = DiscoverySNPCallerPlugin.readTaxaFsFromFile(pedigreeFile);
            if (taxaFs == null) {
                throw new IllegalArgumentException("Problem reading the pedigree file. Progam aborted.");
            }
            usePedigree = true;
        }
        if (myArgsEngine.getBoolean("-mnMAF")) {
            minMAF = Double.parseDouble(myArgsEngine.getString("-mnMAF"));
        }
        if (myArgsEngine.getBoolean("-mxMAF")) {
            maxMAF = Double.parseDouble(myArgsEngine.getString("-mxMAF"));
        }
        if (myArgsEngine.getBoolean("-hLD")) {
            hLD = true;
        }
        if (myArgsEngine.getBoolean("-mnR2")) {
            if (hLD) {
                minR2 = Double.parseDouble(myArgsEngine.getString("-mnR2"));
            } else {
                printUsage();
                throw new IllegalArgumentException("The -mnR2 option requires that the -hLD option is invoked\n");
            }
        }
        if (myArgsEngine.getBoolean("-mnBonP")) {
            if (hLD) {
                minBonP = Double.parseDouble(myArgsEngine.getString("-mnBonP"));
            } else {
                printUsage();
                throw new IllegalArgumentException("The -mnBonP option requires that the -hLD option is invoked\n");
            }
        }
        if (myArgsEngine.getBoolean("-hmp")) {
            if (myArgsEngine.getBoolean("-vcf")){
                throw new IllegalArgumentException("-hmp and -vcf options are mutual exclusive!\n");
            }
            suppliedInputFileName = myArgsEngine.getString("-hmp");
            inputFormat = INPUT_FORMAT.hapmap;
        } else if (myArgsEngine.getBoolean("-vcf")) {
            suppliedInputFileName = myArgsEngine.getString("-vcf");
            inputFormat = INPUT_FORMAT.vcf;
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a HapMap or VCF file to filter.\n");
        }
        if (myArgsEngine.getBoolean("-o")) {
            suppliedOutputFileName = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file name.\n");
        }
        if (myArgsEngine.getBoolean("-maxAlleleVCF")) {
            
            if (! myArgsEngine.getBoolean("-vcf")){
                throw new IllegalArgumentException("-maxAlleleVCF option only works with -vcf input.\n");
            } 
            myMaxNumAlleles = Integer.parseInt(myArgsEngine.getString("-maxAlleleVCF"));
        }
        else
        {
            myMaxNumAlleles = VCFUtil.VCF_DEFAULT_MAX_NUM_ALLELES;
        }
        if (myArgsEngine.getBoolean("-snpLog")) {
            snpLogFileName = myArgsEngine.getString("-snpLog");
        }
        snpLogging = new SNPLogging(snpLogFileName, this.getClass());
    }

    public static String[] getLowCoverageLines(Alignment a, double pCoverage) {
        ArrayList<String> lowLines = new ArrayList<String>();
        for (int i = 0; i < a.getSequenceCount(); i++) {
            int covered = 0;
            for (int j = 0; j < a.getSiteCount(); j++) {
                if (a.getBase(i, j) != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    covered++;
                }
            }
            double propCovered = (double) covered / (double) a.getSiteCount();
            // myLogger.info(a.getTaxaName(i)+":"+propCovered);
            if (propCovered < pCoverage) {
                lowLines.add(a.getIdGroup().getIdentifier(i).getFullName());
            }
        }
        String[] lowL = lowLines.toArray(new String[0]);
        return lowL;
    }

    private String[] getHighExpectedFTaxa(Alignment a) {
        ArrayList<String> highFLines = new ArrayList<String>();
        int nInbredTaxa = 0;
        for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
            String fullTaxonName = a.getFullTaxaName(taxon);
            if (taxaFs.containsKey(fullTaxonName)) {
                if (taxaFs.get(fullTaxonName) >= minF) {
                    highFLines.add(fullTaxonName);
                    nInbredTaxa++;
                }
            }
        }
        myLogger.info(nInbredTaxa + " taxa with an Expected F >= the mnF of " + minF + " were found in the pedigree file (-p option)");
        String[] highF = highFLines.toArray(new String[0]);
        return highF;
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
