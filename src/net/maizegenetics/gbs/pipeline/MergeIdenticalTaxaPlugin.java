package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import javax.swing.ImageIcon;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableVCFAlignment;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.VCFUtil;

import org.apache.log4j.Logger;

/**
 * Basic filters needed for removing bad sites and taxa from GBS pipelines
 * @author edbuckler, Gabriel Rodrigues Alves Margarido
 */
public class MergeIdenticalTaxaPlugin extends AbstractPlugin {

    int startChromosome = 1, endChromosome = 10;
    private ArgsEngine myArgsEngine = null;
    private static final Logger myLogger = Logger.getLogger(MergeIdenticalTaxaPlugin.class);
    private String suppliedInputFileName, suppliedOutputFileName, infile, outfile;
    private double majorityRule = 0.8;
    private boolean makeHetCalls = true;
    String[] lowCoverageTaxa = null;
    
    private static enum INPUT_FORMAT {hapmap, vcf}; //input file format, acceptable values are "hapmap" "vcf" 
    private INPUT_FORMAT inputFormat = INPUT_FORMAT.hapmap;
    private int myMaxNumAlleles;

    public MergeIdenticalTaxaPlugin() {
        super(null, false);
    }

    public MergeIdenticalTaxaPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        for (int chr = startChromosome; chr <= endChromosome; chr++) {
            infile = suppliedInputFileName.replace("+", "" + chr);
            outfile = suppliedOutputFileName.replace("+", "" + chr);
            myLogger.info("Reading: " + infile);
            Alignment a = null;
            
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
            AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, false, true);

            IdGroup idg = a.getIdGroup();
            TreeMap<String, List<String>> sortedIds2 = new TreeMap<String, List<String>>();
            int uniqueTaxa = 0;
            for (int i = 0; i < idg.getIdCount(); i++) {
                List<String> l = sortedIds2.get(idg.getIdentifier(i).getName());
                if (l == null) {
                    sortedIds2.put(idg.getIdentifier(i).getName(), l = new ArrayList<String>());
                    uniqueTaxa++;
                }
                l.add(idg.getIdentifier(i).getFullName());
            }
            IdGroup newGroup = new SimpleIdGroup(uniqueTaxa);
            int index = 0;
            for (List<String> l : sortedIds2.values()) {
                if (l.size() > 1) {
                    newGroup.setIdentifier(index, new Identifier(l.get(0).split(":")[0] + ":MERGE"));
                    System.out.println("To be merged: " + l.size() + ": " + l);
                } else {
                    newGroup.setIdentifier(index, new Identifier(l.get(0)));
                }
                //System.out.println(newGroup.getIdentifier(index).getFullName());
                index++;
            }
            System.out.println("Total taxa:" + idg.getIdCount());
            System.out.println("Unique taxa:" + uniqueTaxa);
            //MutableNucleotideAlignment theMSA = new MutableNucleotideAlignment(newGroup, a.getSiteCount(), a.getLoci());
            MutableNucleotideAlignment theMSA = null;
            if (inputFormat == INPUT_FORMAT.hapmap){
                theMSA = MutableNucleotideAlignment.getInstance(newGroup, a.getSiteCount());
            }
            else if (inputFormat == INPUT_FORMAT.vcf){
                theMSA = MutableVCFAlignment.getInstance(newGroup, a.getSiteCount(),newGroup.getIdCount(), a.getSiteCount(), myMaxNumAlleles);
            }
            for (int s = 0; s < a.getSiteCount(); s++) {
                theMSA.setLocusOfSite(s, a.getLocus(s));
                theMSA.setPositionOfSite(s, a.getPositionInLocus(s));
                if (inputFormat == INPUT_FORMAT.vcf){
                    theMSA.setReferenceAllele(s, a.getReferenceAllele(s));
                    theMSA.setCommonAlleles(s, a.getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE.Depth, s));
                }
                
                //theMSA.setSitePrefix(s, (byte) a.getSNPID(s).charAt(0));
                //theMSA.setStrandOfSite(s, (byte) 1);
            }
            System.out.printf("MergeLine Number Scored PropScored HetSites PropHets%n");
            int newTaxon = -1;
            byte[] calls = null;
            for (Map.Entry<String, List<String>> entry : sortedIds2.entrySet()) {
                if (inputFormat == INPUT_FORMAT.hapmap)
                {
                    if (entry.getValue().size() > 1) {
                        calls = consensusCalls(a, entry.getValue(), makeHetCalls, majorityRule);
                        newTaxon = theMSA.getIdGroup().whichIdNumber(entry.getValue().get(0).split(":")[0] + ":MERGE");
                    } else {
                        int oldTaxon = a.getIdGroup().whichIdNumber(entry.getValue().get(0));
                        calls = a.getBaseRange(oldTaxon, 0, a.getSiteCount() - 1);
                        newTaxon = theMSA.getIdGroup().whichIdNumber(entry.getValue().get(0));
                    }
                    for (int s = 0; s < a.getSiteCount(); s++) {
                        theMSA.setBase(newTaxon, s, calls[s]);
                    }
                    if (entry.getValue().size() > 1) {
                        int known = 0, hets = 0;
                        for (int s = 0; s < a.getSiteCount(); s++) {
                            byte cb = theMSA.getBase(newTaxon, s);
                            if (cb != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                                known++;
                            }
                            if (AlignmentUtils.isHeterozygous(cb)) {
                                hets++;
                            }
                        }
                        double pctK = (double) known / (double) a.getSiteCount();
                        double pctH = (double) hets / (double) a.getSiteCount();
                        System.out.printf("%s %d %d %.3g %d %.3g %n", entry.getKey(), entry.getValue().size(), known, pctK, hets, pctH);
                    }
                }
                else if (inputFormat == INPUT_FORMAT.vcf){
                    if (entry.getValue().size() > 1){
                        newTaxon = theMSA.getIdGroup().whichIdNumber(entry.getValue().get(0).split(":")[0] + ":MERGE");
                        List<String> taxa = entry.getValue();
                        
                        //the return result is a two day array result[x][y]
                        //y: site index
                        //x: x=0: genotype calling; x=1 to max: allele depth
                        short[][] genotypeAndDepth = consensusCallsForVCF(a, taxa, myMaxNumAlleles);
                        for (int s = 0; s < a.getSiteCount(); s++) {
                            theMSA.setBase(newTaxon, s, (byte)genotypeAndDepth[0][s]);
                            short[] mydepth = new short[theMSA.getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE.Depth, s).length];
                            for (int al=0; al<mydepth.length; al++)
                            {
                                mydepth[al] = genotypeAndDepth[al+1][s];
                            }
                            theMSA.setDepthForAlleles(newTaxon, s, mydepth);
                        }

                    }
                    else {
                        int oldTaxon = a.getIdGroup().whichIdNumber(entry.getValue().get(0));                    
                        calls = a.getBaseRange(oldTaxon, 0, a.getSiteCount());
                        newTaxon = theMSA.getIdGroup().whichIdNumber(entry.getValue().get(0));
                        for (int s = 0; s < a.getSiteCount(); s++) {
                            theMSA.setBase(newTaxon, s, calls[s]);

                            theMSA.setDepthForAlleles(newTaxon, s, a.getDepthForAlleles(oldTaxon, s));
                        }
                    }
                        
                }
            }
            if (inputFormat == INPUT_FORMAT.hapmap)
            {
                ExportUtils.writeToHapmap(theMSA, false, outfile, '\t', this);
            }
            else if (inputFormat == INPUT_FORMAT.vcf)
            {
                ExportUtils.writeToVCF(theMSA, outfile, '\t');
            }
        }

        return null;
    }

    public static byte[] consensusCalls(Alignment a, List<String> taxa, boolean callhets, double majority) {
        short[][] siteCnt = new short[2][a.getSiteCount()];
        int[] taxaIndex = new int[taxa.size()];
        for (int t = 0; t < taxaIndex.length; t++) {
            taxaIndex[t] = a.getIdGroup().whichIdNumber(taxa.get(t));
        }
        byte[] calls = new byte[a.getSiteCount()];
        Arrays.fill(calls, Alignment.UNKNOWN_DIPLOID_ALLELE);
        for (int s = 0; s < a.getSiteCount(); s++) {
            byte mjAllele = a.getMajorAllele(s);
            byte mnAllele = a.getMinorAllele(s);
            byte mj=AlignmentUtils.getDiploidValue(mjAllele,mjAllele);
            byte mn=AlignmentUtils.getDiploidValue(mnAllele,mnAllele);
       //     byte[] snpValue = {mj, mn};
            //byte het = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(snpValue);
            byte het = AlignmentUtils.getUnphasedDiploidValue(mjAllele, mnAllele);
            for (int t = 0; t < taxaIndex.length; t++) {
                byte ob = a.getBase(taxaIndex[t], s);
                if (ob == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    continue;
                }
                if (ob == mj) {
                    siteCnt[0][s]++;
                } else if (ob == mn) {
                    siteCnt[1][s]++;
                } else if (AlignmentUtils.isEqual(ob, het)) {
                    siteCnt[0][s]++;
                    siteCnt[1][s]++;
                }
            }
            int totalCnt = siteCnt[0][s] + siteCnt[1][s];
            if (totalCnt == 0) {
                continue;  //no data leave missing
            }
            if ((double) siteCnt[0][s] / (double) totalCnt > majority) {
                calls[s] = mj;
            } else if ((double) siteCnt[1][s] / (double) totalCnt > majority) {
                calls[s] = mn;
            } else if (callhets) {
                calls[s] = het;
            }
        }
        return calls;
    }

    public static short[][] consensusCallsForVCF  (Alignment a, List<String> taxa, int MaxNumAlleles)
    {
        //the return result is a two day array result[x][y]
        //y: site index
        //x: x=0: genotype calling; x=1 to max: allele depth
        short[][] result = new short[MaxNumAlleles+1][a.getSiteCount()];
        for (short[] row: result)
        {
            Arrays.fill(row, (short)0);
        }
        Arrays.fill(result[0], (short)(-1));
        
        int[] taxaIndex = new int[taxa.size()];
        for (int t = 0; t < taxa.size(); t++) {
            taxaIndex[t] = a.getIdGroup().whichIdNumber(taxa.get(t));
        }
        for (int s = 0; s < a.getSiteCount(); s++) {
            byte[] alleles = a.getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE.Depth, s);
            

            int[] alleleDepth = new int[alleles.length];
            Arrays.fill(alleleDepth, 0);

            for (int t = 0; t < taxaIndex.length; t++) {
                short[] myAlleledepth = a.getDepthForAlleles(taxaIndex[t], s);
                for (int al=0; al<myAlleledepth.length; al++)
                {
                    alleleDepth[al] += (int)myAlleledepth[al];
                }
            }
            result[0][s] = VCFUtil.resolveVCFGeno(alleles, alleleDepth);
            
            for (int al=0; al<alleles.length; al++){
                result[al+1][s] = (short)(alleleDepth[al]>32767?32767:alleleDepth[al]);
            }
        }
        return result;
    }
    private void printUsage() {
        myLogger.info(
                "Input format:\n"
                + "-hmp      Input HapMap file; use a plus sign (+) as a wild card character to "
                + "            specify multiple chromosome numbers.\n"
                + "-vcf      Input VCF file. Use a plus sign (+) as a wild card character to specify multiple chromosome numbers. Options -hmp and -vcf are mutual exclusive.\n"
                + "-o        Output HapMap file; use a plus sign (+) as a wild card character to "
                + "            specify multiple chromosome numbers.\n"
                + "-xHet     Exclude heterozygotes calls (default: " + makeHetCalls + ")"
                + "-hetFreq  Cutoff frequency between het vs. homozygote calls (default: " + majorityRule + ")"
                + "-sC       Start chromosome (default 1).\n"
                + "-eC       End chromosome (default 10).\n"
                + "-maxAlleleVCF   Maximum number of alleles allowed in vcf file.\n");
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
            myArgsEngine.add("-xHets", "--excludeHets", false);
            myArgsEngine.add("-hetFreq", "--heterozygoteFreqCutoff", true);
            myArgsEngine.add("-sC", "--startChrom", true);
            myArgsEngine.add("-eC", "--endChrom", true);
            myArgsEngine.add("-maxAlleleVCF", "--maxAlleleVCF", true);
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

        if (myArgsEngine.getBoolean("-hmp")) {
            if (myArgsEngine.getBoolean("-vcf")){
                throw new IllegalArgumentException("-hmp and -vcf options are mutual exclusive!\n");
            }
            suppliedInputFileName = myArgsEngine.getString("-hmp");
            inputFormat = INPUT_FORMAT.hapmap;
        } else if (myArgsEngine.getBoolean("-vcf")) {
            suppliedInputFileName = myArgsEngine.getString("-vcf");
            inputFormat = INPUT_FORMAT.vcf;
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("Please specify a HapMap file or VCF to merge taxa.\n");
        }
        if (myArgsEngine.getBoolean("-o")) {
            suppliedOutputFileName = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file name.\n");
        }
        if (myArgsEngine.getBoolean("-xHets")) {
            makeHetCalls = false;
        }
        if (myArgsEngine.getBoolean("-hetFreq")) {
            majorityRule = Double.parseDouble(myArgsEngine.getString("-hetFreq"));
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
