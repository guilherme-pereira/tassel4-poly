/*
 * MergeDuplicateSNPsPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.File;

import javax.swing.ImageIcon;
import org.apache.log4j.Logger;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.plugindef.AbstractPlugin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MutableVCFAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.util.VCFUtil;

/**
 * This class is intended to be run directly after DiscoverySNPCallerPlugin,
 * using the HapMap file from that step as input.
 *
 * It finds duplicate SNPs in the HapMap file, and merges them if they have the
 * same pair of alleles (not necessarily in the same maj/min order) and if there
 * mismatch rate is no greater than the threshold (-maxMisMat). If -callHets is
 * on, then genotypic disagreements will be called heterozygotes (otherwise set
 * to 'N' = default).
 *
 * By default, any remaining unmerged duplicate SNPs (but not indels) will be
 * deleted. They can be kept by invoking the -kpUnmergDups option.
 *
 * If the germplasm is not fully inbred, and still contains residual
 * heterozygosity (like the maize NAM or IBM populations do) then -callHets
 * should be on and -maxMisMat should be set fairly high (0.1 to 0.2, depending
 * on the amount of heterozygosity).
 *
 * @author jcg233, Gabriel Rodrigues Alves Margarido
 */
public class MergeDuplicateSNPsPlugin extends AbstractPlugin {

    private static Logger myLogger = Logger.getLogger(MergeDuplicateSNPsPlugin.class);
    private static ArgsEngine myArgsEngine = null;
    private String suppliedInputFileName, suppliedOutputFileName, infile, outfile;
    private String snpLogFileName;
    private SNPLogging snpLogging = null;
    private double maxMisMat = 0.05;
    private boolean usePedigree = false;
    private HashMap<String, Double> taxaFs = null;
    private boolean[] useTaxaForCompare = null;
    private int nInbredTaxa = Integer.MIN_VALUE;
    private boolean callHets = false;  // true = when two genotypes disagree at a SNP, call it a heterozygote;  false = set to missing;
    private boolean kpUnmergDups = false;  // keep unmerged SNPs (not indels) in the data file
    private int startChr = 1, endChr = 10;
    private static enum INPUT_FORMAT {hapmap, vcf}; //input file format, acceptable values are "hapmap" "vcf" 
    private INPUT_FORMAT inputFormat = INPUT_FORMAT.hapmap;
    private int myMaxNumAlleles;

    public MergeDuplicateSNPsPlugin() {
        super(null, false);
    }

    public MergeDuplicateSNPsPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
                "\n\nUsage is as follows:\n"
                + "-hmp           Input GBS genotype file (in HapMap format). Use a plus sign (+) as a wild card character to specify multiple chromosome numbers.\n"
                + "-vcf           Input GBS genotype file (in VCF format). Use a plus sign (+) as a wild card character to specify multiple chromosome numbers. Options -hmp and -vcf are mutual exclusive.\n"
                + "-o             Output HapMap file. Use a plus sign (+) as a wild card character to specify multiple chromosome numbers.\n"
                + "-misMat        Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: " + maxMisMat + ")\n"
                + "-p             Pedigree file containing full sample names (or expected names after merging) & expected inbreeding\n"
                + "                 coefficient (F) for each.  Only highly inbred taxa, with F >= 0.8 (e.g., S3 or more), will be used\n"
                + "                 to test if two duplicate SNPs agree with each other (default: use ALL taxa to compare duplicate SNPs)\n"
                + "-callHets      When two genotypes disagree at a SNP, call it a heterozygote (default: " + callHets + " = set to missing)\n"
                + "-kpUnmergDups  When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: " + kpUnmergDups + " = delete them)\n"
                + "-sC             Start chromosome\n"
                + "-eC             End chromosome\n"
                + "-maxAlleleVCF   Maximum number of alleles allowed in vcf file.\n"
                + "-snpLog        SNPs Removed Log file name\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-hmp", "--hmpFile", true);
            myArgsEngine.add("-vcf", "--vcfFile", true);
            myArgsEngine.add("-o", "--outFile", true);
            myArgsEngine.add("-misMat", "--maxMismatchRate", true);
            myArgsEngine.add("-p", "--pedigree-file", true);
            myArgsEngine.add("-callHets", "--callHeterozygotes", false);
            myArgsEngine.add("-kpUnmergDups", "--keepUnmergedDuplicates", false);
            myArgsEngine.add("-sC", "--startChromosome", true);
            myArgsEngine.add("-eC", "--endChromosome", true);
            myArgsEngine.add("-maxAlleleVCF", "--maxAlleleVCF", true);
            myArgsEngine.add("-snpLog", "", true);
        }
        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-hmp")) {
            if (myArgsEngine.getBoolean("-vcf")){
                throw new IllegalArgumentException("-hmp and -vcf options are mutual exclusive!\n");
            } 
            suppliedInputFileName = myArgsEngine.getString("-hmp");
            inputFormat = INPUT_FORMAT.hapmap;
        } else if (myArgsEngine.getBoolean("-vcf")){
            suppliedInputFileName = myArgsEngine.getString("-vcf");
            inputFormat = INPUT_FORMAT.vcf;
        }
        else {
            printUsage();
            throw new IllegalArgumentException("Please specify a HapMap or vcf file to filter.\n");
        }
        if (myArgsEngine.getBoolean("-o")) {
            suppliedOutputFileName = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file name.\n");
        }
        if (myArgsEngine.getBoolean("-misMat")) {
            maxMisMat = Double.parseDouble(myArgsEngine.getString("-misMat"));
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
        if (myArgsEngine.getBoolean("-callHets")) {
            callHets = true;
        }
        if (myArgsEngine.getBoolean("-kpUnmergDups")) {
            kpUnmergDups = true;
        }
        if (myArgsEngine.getBoolean("-sC")) {
            startChr = Integer.parseInt(myArgsEngine.getString("-sC"));
        }
        if (myArgsEngine.getBoolean("-eC")) {
            endChr = Integer.parseInt(myArgsEngine.getString("-eC"));
        }
        if (endChr - startChr < 0) {
            printUsage();
            throw new IllegalArgumentException("Error: The start chromosome is higher than the end chromosome.");
        }
        if (myArgsEngine.getBoolean("-snpLog")) {
            snpLogFileName = myArgsEngine.getString("-snpLog");
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
                
                
        snpLogging = new SNPLogging(snpLogFileName, this.getClass());
    }

    @Override
    public DataSet performFunction(DataSet input) {
        for (int chr = startChr; chr <= endChr; chr++) {
            infile = suppliedInputFileName.replace("+", "" + chr);
            outfile = suppliedOutputFileName.replace("+", "" + chr);
            myLogger.info("Reading: " + infile);
            Alignment a;
            try {
                if (inputFormat == INPUT_FORMAT.hapmap)
                {
                    a = ImportUtils.readFromHapmap(infile, this);
                }
                else if (inputFormat == INPUT_FORMAT.vcf)
                {
                    a = ImportUtils.readFromVCF(infile, this, myMaxNumAlleles);
                }
                else
                {
                     throw new IllegalArgumentException("File format " + inputFormat + " is not recognized!");
                }
                
                myLogger.info("Original Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
                if (usePedigree && !maskNonInbredTaxa(a)) {
                    throw new IllegalArgumentException("Mismatch between taxa names in the input hapmap and pedigree files.");
                }
            } catch (Exception e) {
                myLogger.info("Could not read input hapmap file for chr" + chr + ":\n\t" + infile + "\n\te: " + e + "\n\tSkipping...");
                continue;
            }
            MutableNucleotideAlignment msa = null;
            if (inputFormat == INPUT_FORMAT.hapmap)
            {
                msa = MutableNucleotideAlignment.getInstance(a.getIdGroup(), a.getSiteCount());
            }
            else if (inputFormat == INPUT_FORMAT.vcf)
            {
                 msa = MutableVCFAlignment.getInstance(a.getIdGroup(), a.getSiteCount(), myMaxNumAlleles);
            }
            ArrayList<Integer> samePosAL = new ArrayList<Integer>();
            Integer[] samePos = null;
            int currentPos = a.getPositionInLocus(0);
            for (int s = 0; s < a.getSiteCount(); s++) {  // must be sorted by position, as HapMap files typically are (ImportUtils.readFromHapmap() fails if they aren't)
                int newPos = a.getPositionInLocus(s);
                if (newPos == currentPos) {   // assumes that the strands are all '+' (in DiscoverySNPCallerPlugin(), - strand genos were complemented)
                    samePosAL.add(s);  // collect markers with the same position
                } else {
                    samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
                    if (samePosAL.size() > 1) {  // merge sets of 2 or more markers with the same position and alleles (alleles are not necessarily in the same order, maj/min)
                        if (inputFormat == INPUT_FORMAT.hapmap)
                        {
                            processSNPsWithSamePosition(samePos, a, chr, currentPos, msa);
                        }
                        else if (inputFormat == INPUT_FORMAT.vcf)
                        {
                            processSNPsWithSamePositionForVCF(samePos, a, chr, currentPos, msa);
                        }
                    } else {  // site has a unique position: write its genos to the msa
                        byte[] genos = new byte[a.getSequenceCount()];
                        for (int t = 0; t < a.getSequenceCount(); ++t) {
                            genos[t] = a.getBase(t, samePos[0]);
                        }
                        addSiteToMutableAlignment(chr, currentPos, genos, msa);
                        if (inputFormat==INPUT_FORMAT.vcf)
                        {
                            int lastSiteIndex = msa.getSiteCount() -1;
                            msa.setCommonAlleles(lastSiteIndex, a.getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE.Depth, samePos[0]));
                            msa.setReferenceAllele(lastSiteIndex, a.getReferenceAllele(samePos[0]));
                            for (int tt=0; tt<a.getSequenceCount(); tt++)
                            {
                                msa.setDepthForAlleles(tt, lastSiteIndex, a.getDepthForAlleles(tt, samePos[0]));
                            }
                        }
                    }
                    // start a new collection of markers
                    samePosAL = new ArrayList<Integer>();
                    samePosAL.add(s);
                    currentPos = newPos;
                }
            }
            // finish last site or set of sites
            samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
            if (samePosAL.size() > 1) {
                //processSNPsWithSamePosition(samePos, a, chr, currentPos, msa);
                if (inputFormat==INPUT_FORMAT.hapmap)
                {
                    processSNPsWithSamePosition(samePos, a, chr, currentPos, msa);
                }
                else if (inputFormat==INPUT_FORMAT.vcf)
                {
                    processSNPsWithSamePositionForVCF(samePos, a, chr, currentPos, msa);
                }
            } else {  // site has a unique position: write its genos to the msa
                byte[] genos = new byte[a.getSequenceCount()];
                for (int t = 0; t < a.getSequenceCount(); ++t) {
                    genos[t] = a.getBase(t, samePos[0]);
                }
                addSiteToMutableAlignment(chr, currentPos, genos, msa);
                if (inputFormat==INPUT_FORMAT.vcf)
                {
                    int lastSiteIndex = msa.getSiteCount() -1;
                    msa.setCommonAlleles(lastSiteIndex, a.getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE.Depth, samePos[0]));
                    msa.setReferenceAllele(lastSiteIndex, a.getReferenceAllele(samePos[0]));
                    for (int tt=0; tt<a.getSequenceCount(); tt++)
                    {
                        msa.setDepthForAlleles(tt, lastSiteIndex, a.getDepthForAlleles(tt, samePos[0]));
                    }
                }
            }
            msa.clean();
            myLogger.info("Number of sites written after merging duplicate SNPs: " + msa.getSiteCount());
            
            if (inputFormat==INPUT_FORMAT.hapmap)
            {
                if (!kpUnmergDups) {
                    deleteRemainingDuplicates(msa);
                }

                ExportUtils.writeToHapmap(msa, false, outfile, '\t', this);
            }
            else if (inputFormat==INPUT_FORMAT.vcf)
            {
                ExportUtils.writeToVCF(msa, outfile, '\t');
            }
        }
        return null;
    }

    private void processSNPsWithSamePosition(Integer[] samePos, Alignment a, int chr, int currentPos, MutableNucleotideAlignment msa) {
        boolean[] finished = new boolean[samePos.length];
        for (int i = 0; i < finished.length; ++i) {
            finished[i] = false;   // indicates if the site has already been merged with a previous site OR written as is to msa
        }
        for (int s1 = 0; s1 < samePos.length - 1; ++s1) {
            if (finished[s1]) {
                continue;  // s1 has already been merged with a previous site
            }
            byte[] currMerge = new byte[a.getSequenceCount()];
            for (int t = 0; t < a.getSequenceCount(); ++t) {
                currMerge[t] = a.getBase(t, samePos[s1]); // set the current merger of genotypes to those for site s1
            }
            byte[] currAlleles = new byte[2];
            currAlleles[0] = a.getMajorAllele(samePos[s1].intValue());
            currAlleles[1] = a.getMinorAllele(samePos[s1].intValue());
            if (currAlleles[0] == NucleotideAlignmentConstants.GAP_ALLELE || currAlleles[1] == NucleotideAlignmentConstants.GAP_ALLELE) {
                addSiteToMutableAlignment(chr, currentPos, currMerge, msa);
                finished[s1] = true;
                continue;
            }
            Arrays.sort(currAlleles);
            for (int s2 = s1 + 1; s2 < samePos.length; ++s2) {
                if (finished[s2]) {
                    continue;  // s2 has already been merged with a previous site (perhaps with different alleles)
                }
                byte[] newAlleles = new byte[2];
                newAlleles[0] = a.getMajorAllele(samePos[s2].intValue());
                newAlleles[1] = a.getMinorAllele(samePos[s2].intValue());
                Arrays.sort(newAlleles);
                if (newAlleles[0] == NucleotideAlignmentConstants.GAP_ALLELE || newAlleles[1] == NucleotideAlignmentConstants.GAP_ALLELE) {
                    continue;
                }
                // Check if the alleles match.  If they do, merge the genos, provided that the number of genotypic mismatches is below threshold
                if (Arrays.equals(currAlleles, newAlleles)) {
                    int nMismatch = 0;
                    int nCompare = 0;
                    byte[] possibleMerge = new byte[a.getSequenceCount()];
                    for (int t = 0; t < a.getSequenceCount(); ++t) {
                        byte geno2 = a.getBase(t, samePos[s2]);
                        if (currMerge[t] != Alignment.UNKNOWN_DIPLOID_ALLELE && geno2 != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            if (!usePedigree || useTaxaForCompare[t]) {
                                ++nCompare;
                            }
                            if (!AlignmentUtils.isEqual(currMerge[t], geno2)) {
                                if (!usePedigree || useTaxaForCompare[t]) {
                                    ++nMismatch;
                                }
                                try {
                                    possibleMerge[t] = callHets ? resolveHet(currMerge[t], geno2) : Alignment.UNKNOWN_DIPLOID_ALLELE;
                                } catch (Exception e) {
                                    myLogger.warn(
                                            "Invalid genotypes (" + a.getBaseAsString(t, samePos[s1]) + " and " + a.getBaseAsString(t, samePos[s2]) + ") at position:" + currentPos + " taxon:" + a.getTaxaName(t));
                                }
                            } else {
                                possibleMerge[t] = currMerge[t];
                            }
                        } else if (currMerge[t] == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            possibleMerge[t] = geno2;
                        } else if (geno2 == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            possibleMerge[t] = currMerge[t];
                        }
                    }
                    if ((nCompare == 0) || ((double) nMismatch / nCompare <= maxMisMat)) {
                        for (int t = 0; t < a.getSequenceCount(); ++t) {
                            currMerge[t] = possibleMerge[t];
                        }
                        double valueMisMat = Double.NaN;
                        if (nCompare != 0) {
                            valueMisMat = (double) nMismatch / nCompare;
                        }
                        snpLogging.writeEntry(a, samePos[s2], null, null, "Genotypes Less Max Mismatch", "Merged", String.valueOf(valueMisMat), String.valueOf(maxMisMat));
                        finished[s2] = true;
                    }
                }
            }
            // Finished comparing s1 to all other sites
            // Write the current "merged" genos to the msa (if all comparisons disagreed, these will be the original s1 genos)
            addSiteToMutableAlignment(chr, currentPos, currMerge, msa);
            finished[s1] = true;
        }
        // make sure that everything in the current collection (samePos) got finished
        //  (for example, the final site may not have been merged with anything, or might have had gap as an allele)
        for (int site = 0; site < finished.length; ++site) {
            if (!finished[site]) {
                byte[] genos = new byte[a.getSequenceCount()];
                for (int t = 0; t < a.getSequenceCount(); ++t) {
                    genos[t] = a.getBase(t, samePos[site]);
                }
                addSiteToMutableAlignment(chr, currentPos, genos, msa);
            }
        }
    }

   private void processSNPsWithSamePositionForVCF(Integer[] samePos, Alignment a, int chr, int currentPos, MutableNucleotideAlignment msa)
    {
        int taxaCount = a.getSequenceCount();
        //merged depth table
        HashMap <Byte,int[]> MergedDepthTable = new HashMap<Byte,int[]>();
        HashMap <Byte,Integer> AllAlleleDepthTable = new HashMap<Byte,Integer>();
        
        for (int s:samePos)
        {
            //get alleles for each site
            byte[] alleles = a.getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE.Depth, s);
            
            //initiate new alleles in the allele list into the hashmap
            for (byte allele:alleles)
            {
                if (!MergedDepthTable.containsKey(allele))
                {
                    int[] newDepth = new int[taxaCount];
                    Arrays.fill(newDepth, 0);
                    MergedDepthTable.put(allele, newDepth);
                    AllAlleleDepthTable.put(allele, 0);
                }             
            }
            //add up the depth for each alleles
            for (int t = 0; t < a.getSequenceCount(); ++t) 
            {
                short[] currDepth = a.getDepthForAlleles(t, s); // get the current depth for site s taxa t 
                for (int i=0; i < alleles.length; i++)
                {
                    if (currDepth[i]>0)
                    {
                        MergedDepthTable.get(alleles[i])[t] +=  currDepth[i];
                        AllAlleleDepthTable.put(alleles[i], AllAlleleDepthTable.get(alleles[i]) + currDepth[i]);
                    }
                }
            }
        }           
        
        //sort allele by count
        ArrayList <Byte> AllAllelesList = new ArrayList <Byte>();
        AllAllelesList.addAll(AllAlleleDepthTable.keySet());
        Collections.sort(AllAllelesList, new HashValueComparator(AllAlleleDepthTable));
        
        //get top alleles and read depth for each allele
        int allelesCount = (AllAllelesList.size()>myMaxNumAlleles)?myMaxNumAlleles:AllAllelesList.size();
        byte CommonAlleles[] = new byte[allelesCount]; //all alleles
        int[][] alleleDepthsInTaxa = new int[myMaxNumAlleles][taxaCount];
        for(int[] row: alleleDepthsInTaxa)
        {
            Arrays.fill(row, 0);
        }
        
        for (int i=0; i<allelesCount; i++)
        {
            CommonAlleles[i] = AllAllelesList.get(i);
            for (int t=0; t<taxaCount; t++)
            {
                alleleDepthsInTaxa[i][t] = MergedDepthTable.get(CommonAlleles[i])[t];
            }
        }
        
        //resolve genotypies
        byte[] genos = new byte[taxaCount];
        for (int t=0; t<taxaCount; t++)
        {
               genos[t] = VCFUtil.resolveVCFGeno(CommonAlleles,alleleDepthsInTaxa, t);
        }
        
        //calculate mismatch rate
        int nCompared =0;
        int nMisMatch =0;
        for (int t=0; t<genos.length; t++)
        {
            byte mergedGeno = genos[t];
            if (mergedGeno == Alignment.UNKNOWN_DIPLOID_ALLELE)
            {
                continue;
            }
            for (int s:samePos)
            {
                byte singleBase = a.getBase(t, s);
                if (singleBase==Alignment.UNKNOWN_DIPLOID_ALLELE)
                {
                    continue;
                }
                nCompared ++;
                if (!AlignmentUtils.isEqual(mergedGeno, singleBase))
                {
                    nMisMatch ++;
                }
            }
        }
        double myMisMatchRate=0;
        if (nCompared>0)
        {
            myMisMatchRate=(double)nMisMatch/nCompared;
        }
        if (myMisMatchRate<maxMisMat)
        {
            addSiteToMutableAlignment(chr, currentPos, genos, msa);
            int lastSiteIndex = msa.getSiteCount() - 1;
            msa.setCommonAlleles(lastSiteIndex, CommonAlleles);
            msa.setReferenceAllele(lastSiteIndex, a.getReferenceAllele(samePos[0]));
            for (int t=0; t<taxaCount; t++)
            {
                short[] alleleDepth = new short[allelesCount];
                for (int i=0; i<allelesCount; i++)
                {
                    alleleDepth[i] = alleleDepthsInTaxa[i][t]>32767?(short)32767:(short)alleleDepthsInTaxa[i][t];
                
                }
                msa.setDepthForAlleles(t, lastSiteIndex, alleleDepth);
            }
        }
        else 
        {
            System.out.println("Not merged position: " + a.getPositionInLocus(samePos[0]) +  " Mismatch: "+ (int)(myMisMatchRate *100) + "%." );
            
            
            if (kpUnmergDups)
            {
                for (int s:samePos)
                {
                    genos = new byte[a.getSequenceCount()];
                    for (int t = 0; t < a.getSequenceCount(); ++t) {
                        genos[t] = a.getBase(t, s);
                    }
                    addSiteToMutableAlignment(chr, currentPos, genos, msa);

                    int lastSiteIndex = msa.getSiteCount() -1;
                    msa.setCommonAlleles(lastSiteIndex, a.getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE.Depth, s));
                    msa.setReferenceAllele(lastSiteIndex, a.getReferenceAllele(s));
                    for (int tt=0; tt<a.getSequenceCount(); tt++)
                    {
                        msa.setDepthForAlleles(tt, lastSiteIndex, a.getDepthForAlleles(tt, s));
                    }

                }
            }
        }
    }
    
    private void addSiteToMutableAlignment(int chromosome, int position, byte[] genos, MutableNucleotideAlignment theMSA) {
        int currSite = theMSA.getSiteCount();
        //int currSite = theMSA.getNextFreeSite();
        theMSA.addSite(currSite);
        theMSA.setLocusOfSite(currSite, new Locus(String.valueOf(chromosome), String.valueOf(chromosome), -1, -1, null, null));
        theMSA.setPositionOfSite(currSite, position);
        //theMSA.setStrandOfSite(currSite, (byte) '+');
        for (int tx = 0; tx < genos.length; tx++) {
            theMSA.setBase(tx, currSite, genos[tx]);
        }
    }

    private void deleteRemainingDuplicates(MutableNucleotideAlignment theMSA) {

        ArrayList<Integer> samePosAL = new ArrayList<Integer>();
        int currentPos = theMSA.getPositionInLocus(0);
        for (int s = 0; s < theMSA.getSiteCount(); s++) {
            int newPos = theMSA.getPositionInLocus(s);
            if (newPos == currentPos) {
                samePosAL.add(s);  // collect markers with the same position
            } else {
                if (samePosAL.size() > 1) {
                    Integer[] samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
                    for (int i = 0; i < samePos.length; ++i) {
                        if (theMSA.getMajorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE
                                && theMSA.getMinorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {
                            snpLogging.writeEntry(theMSA, samePos[i], null, null, "Delete Remaining Duplicates", "Removed", null, null);
                            theMSA.clearSiteForRemoval(samePos[i]);
                        }
                    }
                }
                // start a new collection of markers
                samePosAL = new ArrayList<Integer>();
                samePosAL.add(s);
                currentPos = newPos;
            }
        }
        // finish last site or set of sites
        if (samePosAL.size() > 1) {
            Integer[] samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
            for (int i = 0; i < samePos.length; ++i) {
                if (theMSA.getMajorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE
                        && theMSA.getMinorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {
                    snpLogging.writeEntry(theMSA, samePos[i], null, null, "Delete Remaining Duplicates", "Removed", null, null);
                    theMSA.clearSiteForRemoval(samePos[i]);
                }
            }
        }
        snpLogging.close();
        theMSA.clean();
        myLogger.info("Number of sites written after deleting any remaining, unmerged duplicate SNPs: " + theMSA.getSiteCount());
    }

    private boolean maskNonInbredTaxa(Alignment a) {
        useTaxaForCompare = new boolean[a.getSequenceCount()];  // initialized to false
        nInbredTaxa = 0;
        try {
            for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
                if (taxaFs.containsKey(a.getFullTaxaName(taxon))) {
                    if (taxaFs.get(a.getFullTaxaName(taxon)) >= 0.8) {
                        useTaxaForCompare[taxon] = true;
                        nInbredTaxa++;
                    }
                } else {
                    if (a.getFullTaxaName(taxon).contentEquals("REFERENCE_GENOME")) {
                        useTaxaForCompare[taxon] = false;
                    } else {
                        throw new Exception("Taxon " + a.getFullTaxaName(taxon) + " not found in the pedigree file");
                    }
                }
            }
            myLogger.info(nInbredTaxa + " highly inbred taxa (with an expected F >= 0.8) were found in the input hapmap file (according to the pedigree file)");
            return true;
        } catch (Exception e) {
            myLogger.error("Mismatch between taxa names in the input hapmap file and the pedigree file e=" + e);
            e.printStackTrace();
            return false;
        }
    }

    private static byte resolveHet(byte geno1, byte geno2) {

        byte[] result = new byte[2];
        result[0] = (byte) (geno1 >>> 4);
        byte temp = (byte) (geno1 & 0xf);
        int count = 1;
        if (temp != result[0]) {
            result[count++] = temp;
        }

        temp = (byte) (geno2 >>> 4);
        if (temp == result[0]) {
            // do nothing
        } else if (count == 1) {
            result[count++] = temp;
        } else if (temp != result[1]) {
            throw new IllegalStateException();
        }

        temp = (byte) (geno2 & 0xf);
        if (temp == result[0]) {
            // do nothing
        } else if (count == 1) {
            result[count++] = temp;
        } else if (temp != result[1]) {
            throw new IllegalStateException();
        }

        return (byte) ((result[0] << 4) | (result[1]));

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
class HashValueComparator implements Comparator {

  Map base;
  public HashValueComparator(Map base) {
      this.base = base;
  }

  public int compare(Object a, Object b) {

    if((Integer)base.get(a) < (Integer)base.get(b)) {
      return 1;
    } else if((Integer)base.get(a) == (Integer)base.get(b)) {
      return 0;
    } else {
      return -1;
    }
  }
}