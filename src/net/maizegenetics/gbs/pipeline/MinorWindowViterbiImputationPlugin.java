/*
 * MinorWindowViterbiImputationPlugin
 */
package net.maizegenetics.gbs.pipeline;


import java.awt.Frame;
import java.io.File;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.ImageIcon;
import net.maizegenetics.annotation.Citation;
import net.maizegenetics.gwas.imputation.EmissionProbability;
import net.maizegenetics.gwas.imputation.TransitionProbability;
import net.maizegenetics.gwas.imputation.ViterbiAlgorithm;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.alignment.ProjectionAlignment;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.DonorHypoth;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;

import org.apache.log4j.Logger;

/**
 * Imputation approach that relies on nearest neighbor searches of defined haplotypes, 
 * followed by HMM Viterbi resolution or block-based resolution.
 * In June 2013, it is the best approach for substantially unrelated taxa.
 * <p>
 * The algorithm relies on having a series of donor haplotypes.  If phased haplotypes are 
 * already known, they can be used directly.  If they need to be reconstructed, the 
 * {@link FindMergeHaplotypesPlugin} can be used to create haplotypes for windows
 * across the genome.
 * <p>
 * Imputation is done one taxon at a time using the donor haplotypes.  The strategy is as
 *  follows:
 * <li> Every 64 sites is considered a block (or a word in the bitset terminology - length of a long).  There is
 * a separate initial search for nearest neighbors centered around each block. The flanks
 * of the scan window vary but always contain the number of minor alleles specified (a key parameter).
 * Minor alleles approximate the information content of the search.
 * <li> Calculate distance between the target taxon and the donor haplotypes for every window, and rank
 * the top 10 donor haplotypes for each 64 site focus block.
 * <li> Evaluate by Viterbi whether 1 or 2 haplotypes will explain all the sites across the 
 * donor haplotype window.  If successful, move to the next region for donor haplotypes.
 * <li> Resolve each focus block by nearest neighbors.  If inbred NN are not close enough,
 * then do a hybrid search seeded with one parent being the initial 10 haplotypes.  
 * Set bases based on what is resolved first inbred or hybrid.
 *  
 *<p>
 * Error rates are bounded away from zero, but adding 0.5 error to all error
 * rates that that were observed to be zero.
 * <p>
 * TODO:
 * <li>Smaller scale Viterbi, evaluate whether there are regions where 1 or 2 haplotypes are very common (>50% blocks)
 * <li>Use a diversity index of the DHs to determine what transition/emmission matrices to use
 * <li>Move accuracy to one method outside of setAlignmentWithDonors
 *
 * @author Edward Buckler
 * @cite 
 */
@Citation("Two papers: Viterbi from Bradbury, et al (in prep) Recombination patterns in maize\n"+
        "NearestNeighborSearch Swarts,...,Buckler (in prep) Imputation with large genotyping by sequencing data\n")
public class MinorWindowViterbiImputationPlugin extends AbstractPlugin {
    private int startChr, endChr;
    private String hmpFile;
    private String donorFile;
    private String outFileBase;
    private String errFile=null;
    private boolean inbredNN=true;
    private boolean hybridNN=true;
    private int minMinorCnt=20;
    private int minMajorRatioToMinorCnt=10;  //refinement of minMinorCnt to account for regions with all major
    private int maxDonorHypotheses=10;  //number of hypotheses of record from an inbred or hybrid search of a focus block
    private boolean isOutputProjection=false;
    private boolean isChromosomeMode=false;  //deprecated approach to run one chromosome at a time to separate files
    
    private double maximumInbredError=0.02;  //inbreds are tested first, if too much error hybrids are tested.
    private double maxHybridErrorRate=0.005;
    private int minTestSites=100;  //minimum number of compared sites to find a hit
    
    
    private Alignment unimpAlign;  //the unimputed alignment to be imputed, unphased
    private int testing=0;  //level of reporting to stdout
    //major and minor alleles can be differ between the donor and unimp alignment 
    private boolean isSwapMajorMinor=true;  //if swapped try to fix it
    
    private boolean resolveHetIfUndercalled=true;
    //variables for tracking accuracy 
    private int totalRight=0, totalWrong=0, totalHets=0; //global variables tracking errors on the fly
    private int[] siteErrors, siteCorrectCnt, taxonErrors, taxonCorrectCnt;  //error recorded by sites
    
    
    //initialize the transition matrix
    double[][] transition = new double[][] {
                    {.999,.0001,.0003,.0001,.0005},
                    {.0002,.999,.00005,.00005,.0002},
                    {.0002,.00005,.999,.00005,.0002},
                    {.0002,.00005,.00005,.999,.0002},
                    {.0005,.0001,.0003,.0001,.999}
        };
    //initialize the emission matrix, states (5) in rows, observations (3) in columns
    double[][] emission = new double[][] {
                    {.98,.001,.001},
                    {.6,.2,.2},
                    {.4,.2,.4},
                    {.2,.2,.6},
                    {.001,.001,.98}
    };
    
    
//    int[] donorIndices;
    private static ArgsEngine engine = new ArgsEngine();
    private static final Logger myLogger = Logger.getLogger(MinorWindowViterbiImputationPlugin.class);


    public MinorWindowViterbiImputationPlugin() {
        super(null, false);
    }

    public MinorWindowViterbiImputationPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    /**
     * 
     * @param donorFile should be phased haplotypes 
     * @param unImpTargetFile sites must match exactly with donor file
     * @param exportFile output file of imputed sites
     * @param minMinorCnt determines the size of the search window, low recombination 20-30, high recombination 10-15
     * @param hybridMode true=do hybrid search
     */
    public void runMinorWindowViterbiImputation(String donorFile, String unImpTargetFile, 
            String exportFile, int minMinorCnt, int minTestSites, int minSitesPresent, 
            double maxHybridErrorRate, boolean isOutputProjection, boolean imputeDonorFile) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        System.out.println("Retain Rare alleles is:"+TasselPrefs.getAlignmentRetainRareAlleles());
        Alignment[] donorAlign=loadDonors(donorFile);
        this.minTestSites=minTestSites;
        this.isOutputProjection=isOutputProjection;
        unimpAlign=ImportUtils.readGuessFormat(unImpTargetFile, false);
//        if(unImpTargetFile.contains(".h5")) {
//            unimpAlign=BitAlignmentHDF5.getInstance(unImpTargetFile, false);
//        } else {
//            unimpAlign=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
//        }
        unimpAlign.optimizeForTaxa(null);
        
        OpenBitSet[][] conflictMasks=createMaskForAlignmentConflicts(unimpAlign, donorAlign, true);   

        siteErrors=new int[unimpAlign.getSiteCount()];
        siteCorrectCnt=new int[unimpAlign.getSiteCount()];
        taxonErrors=new int[unimpAlign.getSequenceCount()];
        taxonCorrectCnt=new int[unimpAlign.getSequenceCount()];

        System.out.printf("Unimputed taxa:%d sites:%d %n",unimpAlign.getSequenceCount(),unimpAlign.getSiteCount());       
        System.out.println("Creating mutable alignment");
        MutableAlignment mna=null;
        if(isOutputProjection) {
            mna=new ProjectionAlignment(donorAlign[0], unimpAlign.getIdGroup());
        } else {
            if(exportFile.contains("hmp.h5")) {
                ExportUtils.writeToMutableHDF5(unimpAlign, exportFile, new SimpleIdGroup(0), false);
                mna=MutableNucleotideAlignmentHDF5.getInstance(exportFile);
            }else {
                mna=MutableNucleotideAlignment.getInstance(this.unimpAlign);
            }

        }
        long time=System.currentTimeMillis();
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService pool = Executors.newFixedThreadPool(numThreads);
        for (int taxon = 0; taxon < unimpAlign.getSequenceCount(); taxon+=1) {
            ImputeOneTaxon theTaxon=new ImputeOneTaxon(taxon, donorAlign, minSitesPresent, conflictMasks,imputeDonorFile, mna);
    //        theTaxon.run();
            pool.execute(theTaxon);
        }
        pool.shutdown();
        try{
            if (!pool.awaitTermination(48, TimeUnit.HOURS)) {
                System.out.println("processing threads timed out.");
            }
        }catch(Exception e) {
            System.out.println("Error processing threads");
        }
        System.out.println("");
        System.out.println("Time:"+(time-System.currentTimeMillis()));
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile, unImpTargetFile, minMinorCnt));
        System.out.println(s.toString());
        double errRate=(double)totalWrong/(double)(totalRight+totalWrong);
        System.out.printf("TotalRight %d  TotalWrong %d TotalHets: %d ErrRateExcHet:%g %n",totalRight, totalWrong, totalHets, errRate);
//        for (int i = 0; i < siteErrors.length; i++) {
//            System.out.printf("%d %d %d %g %g %n",i,siteCallCnt[i],siteErrors[i], 
//                    (double)siteErrors[i]/(double)siteCallCnt[i], unimpAlign.getMinorAlleleFrequency(i));
//        }
        if(isOutputProjection) {
           ((ProjectionAlignment)mna).save(exportFile);
        } else {
            if(mna instanceof MutableNucleotideAlignmentHDF5) {
                mna.clean();
            } else {
                ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);}
        }
        System.out.printf("%d %g %d %n",minMinorCnt, maximumInbredError, maxDonorHypotheses);
        
    }
    
    private class ImputeOneTaxon implements Runnable{
        int taxon;
        Alignment[] donorAlign;
        int minSitesPresent;
        OpenBitSet[][] conflictMasks;
        boolean imputeDonorFile;
        MutableAlignment mna;
        
        
        public ImputeOneTaxon(int taxon, Alignment[] donorAlign, int minSitesPresent, OpenBitSet[][] conflictMasks,
            boolean imputeDonorFile, MutableAlignment mna) {
            this.taxon=taxon;
            this.donorAlign=donorAlign;
            this.minSitesPresent=minSitesPresent;
            this.conflictMasks=conflictMasks;
            this.imputeDonorFile=imputeDonorFile;
            this.mna=mna;
        }

        
        
        @Override
        public void run() {
            StringBuilder sb=new StringBuilder();
            String name=unimpAlign.getIdGroup().getIdentifier(taxon).getFullName();
            ImputedTaxon impTaxon=new ImputedTaxon(taxon, unimpAlign.getBaseRow(taxon),isOutputProjection);
            int[] unkHets=countUnknownAndHets(impTaxon.getOrigGeno());
            sb.append(String.format("Imputing %d:%s Mj:%d, Mn:%d Unk:%d Hets:%d... ", taxon,name,
                    unimpAlign.getAllelePresenceForAllSites(taxon, 0).cardinality(),
                    unimpAlign.getAllelePresenceForAllSites(taxon, 1).cardinality(), unkHets[0], unkHets[1]));
            boolean enoughData=(unimpAlign.getTotalNotMissingForTaxon(taxon)>minSitesPresent);
//                System.out.println("Too much missing data");
//                continue;
//            }
            int countFullLength=0;
            for (int da = 0; (da < donorAlign.length)&&enoughData ; da++) {
                int donorOffset=unimpAlign.getSiteOfPhysicalPosition(donorAlign[da].getPositionInLocus(0), donorAlign[da].getLocus(0));
                int blocks=donorAlign[da].getAllelePresenceForAllSites(0, 0).getNumWords();
                BitSet[] maskedTargetBits=arrangeMajorMinorBtwAlignments(unimpAlign, taxon, donorOffset, 
                        donorAlign[da].getSiteCount(),conflictMasks[da][0],conflictMasks[da][1]); 
                int[] donorIndices;
                if(imputeDonorFile){
                    donorIndices=new int[donorAlign[da].getSequenceCount()-1];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i; if(i>=taxon) donorIndices[i]++;}
                } else {
                    donorIndices=new int[donorAlign[da].getSequenceCount()];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i;}
                }
                DonorHypoth[][] regionHypth=new DonorHypoth[blocks][maxDonorHypotheses];
                calcInbredDist(impTaxon,maskedTargetBits, donorAlign[da]);
                for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                    int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),maskedTargetBits[1].getBits(), focusBlock, minMinorCnt);
                    if(resultRange==null) continue; //no data in the focus Block
                    //search for the best inbred donors for a segment
                    regionHypth[focusBlock]=getBestInbredDonors(taxon, impTaxon, resultRange[0],resultRange[2], focusBlock, donorAlign[da], donorIndices);
                }
                impTaxon.setSegmentSolved(false);
                //tries to solve the entire donorAlign region with 1 or 2 donor haplotypes
                //consider also running it backwards
                impTaxon=apply1or2Haplotypes(taxon, donorAlign[da], donorOffset, regionHypth,  impTaxon, maskedTargetBits, maxHybridErrorRate);
                if(impTaxon.isSegmentSolved()) {
//                    System.out.printf("VertSolved da:%d L:%s%n",da, donorAlign[da].getLocus(0));
                    countFullLength++; continue;}
                //Kelly try to do something better
           //     impTaxon=applyKellyHaplotypes();
                if(impTaxon.isSegmentSolved()) {countFullLength++; continue;}
                //resorts to solving block by block, first by inbred, and then by hybrid
                if(inbredNN) {
                    impTaxon=applyBlockNN(impTaxon, taxon, donorAlign[da], donorOffset, regionHypth, hybridNN, maskedTargetBits, 
                            minMinorCnt, maxHybridErrorRate, donorIndices);
                }   
            }          
            sb.append(String.format("Viterbi:%d BlocksSolved:%d ", countFullLength, impTaxon.getBlocksSolved()));
            int[] unk=countUnknownAndHets(impTaxon.resolveGeno);
            sb.append(String.format("Unk:%d PropMissing:%g ", unk[0], (double)unk[0]/(double)mna.getSiteCount()));
            sb.append(String.format("Het:%d PropHet:%g ", unk[1], (double)unk[1]/(double)mna.getSiteCount()));
            if(!isOutputProjection) {
                if(mna instanceof MutableNucleotideAlignmentHDF5) {
                    Identifier tID=unimpAlign.getIdGroup().getIdentifier(taxon);
//                    mna.addTaxon(tID);
//                    int newTaxaIndex=mna.getIdGroup().whichIdNumber(tID);
// //                   if(taxon!=newTaxaIndex) System.out.println(taxon+":"+newTaxaIndex);
//                    ((MutableNucleotideAlignmentHDF5)mna).setAllBases(newTaxaIndex, impTaxon.resolveGeno);
                    ((MutableNucleotideAlignmentHDF5)mna).addTaxon(tID, impTaxon.resolveGeno, null);
                } else {
                    mna.setBaseRange(taxon, 0, impTaxon.resolveGeno);
                }
                    
            }
            double errRate=calcErrorForTaxonAndSite(impTaxon); 
            sb.append(String.format("ErR:%g ", errRate));
//            double rate=(double)taxon/(double)(System.currentTimeMillis()-time);
//            double remaining=(unimpAlign.getSequenceCount()-taxon)/(rate*1000);
//            System.out.printf("TimeLeft:%.1fs %n", remaining);
            System.out.println(sb.toString());
            if(isOutputProjection) {
               // System.out.println(breakPoints.toString());
                ((ProjectionAlignment)mna).setCompositionOfTaxon(taxon, impTaxon.breakPoints);
            }
        }
        
    }
    
     public Alignment[] loadDonors(String donorFileRoot){
        File theDF=new File(donorFileRoot);
        String prefilter=theDF.getName().split(".gX.")[0]+".gc"; //grabs the left side of the file
        String prefilterOld=theDF.getName().split("s\\+")[0]+"s"; //grabs the left side of the file
        ArrayList<File> d=new ArrayList<File>();
        for (File file : theDF.getParentFile().listFiles()) {
            if(file.getName().equals(theDF.getName())) {d.add(file);}
             if(file.getName().startsWith(prefilter)) {d.add(file);}
             if(file.getName().startsWith(prefilterOld)) {d.add(file);}
         }
        Alignment[] donorAlign=new Alignment[d.size()];
        for (int i = 0; i < donorAlign.length; i++) {
            System.out.println("Starting Read");
            donorAlign[i]=ImportUtils.readFromHapmap(d.get(i).getPath(), true, (ProgressListener)null);     
            System.out.printf("Donor file:%s taxa:%d sites:%d %n",d.get(i).getPath(), donorAlign[i].getSequenceCount(),donorAlign[i].getSiteCount());
            donorAlign[i].optimizeForTaxa(null);
            System.out.println("Taxa Optimization Done");
            //createMaskForAlignmentConflicts(unimpAlign,donorAlign[i],true);
        }
    	return donorAlign;
    }

    private ImputedTaxon apply1or2Haplotypes(int taxon, Alignment donorAlign, int donorOffset, 
            DonorHypoth[][] regionHypth, ImputedTaxon impT,
            BitSet[] maskedTargetBits, double maxHybridErrorRate) {
        
        int blocks=maskedTargetBits[0].getNumWords();
        //do flanking search 
        if(testing==1) System.out.println("Starting complete hybrid search");
        int[] d=getAllBestDonorsAcrossChromosome(regionHypth,blocks/20);  //TODO 
        DonorHypoth[] best2donors=getBestHybridDonors(taxon, maskedTargetBits[0].getBits(),
                maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, donorAlign, d, d, true);
        if(testing==1) System.out.println(Arrays.toString(best2donors));
        ArrayList<DonorHypoth> goodDH=new ArrayList<DonorHypoth>();
        for (DonorHypoth dh : best2donors) {
            if((dh!=null)&&(dh.getErrorRate()<maxHybridErrorRate)) {
                if(dh.isInbred()==false){
                    dh=getStateBasedOnViterbi(dh, donorOffset, donorAlign);
                }
                if(dh!=null) goodDH.add(dh);
            }
        }
        if(goodDH.isEmpty()) return impT;
        DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
        for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
        impT.setSegmentSolved(true);
        return setAlignmentWithDonors(donorAlign,vdh, donorOffset,false,impT);
    }
    
    private ImputedTaxon applyKellyHaplotypes(int taxon, Alignment donorAlign, int donorOffset, 
            DonorHypoth[][] regionHypth, ImputedTaxon impT,
            BitSet[] maskedTargetBits, double maxHybridErrorRate) {
        
        int blocks=maskedTargetBits[0].getNumWords();
        //do flanking search 
        if(testing==1) System.out.println("Starting complete hybrid search");
        int[] d=getAllBestDonorsAcrossChromosome(regionHypth,5);
        DonorHypoth[] best2donors=getBestHybridDonors(taxon, maskedTargetBits[0].getBits(),
                maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, donorAlign, d, d, true);
        if(testing==1) System.out.println(Arrays.toString(best2donors));
        ArrayList<DonorHypoth> goodDH=new ArrayList<DonorHypoth>();
        for (DonorHypoth dh : best2donors) {
            if((dh!=null)&&(dh.getErrorRate()<maxHybridErrorRate)) {
                if(dh.isInbred()==false){
                    dh=getStateBasedOnViterbi(dh, donorOffset, donorAlign);
                }
                if(dh!=null) goodDH.add(dh);
            }
        }
        if(goodDH.isEmpty()) return impT;
        DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
        for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
        impT.setSegmentSolved(true);
        return setAlignmentWithDonors(donorAlign,vdh, donorOffset,false,impT);
    }  
        
    /**
     * Create mask for all sites where major & minor are swapped in alignments
     * returns in this order [goodMask, swapMjMnMask, errorMask, invariantMask]
     */
    private OpenBitSet[][] createMaskForAlignmentConflicts(Alignment unimpAlign, Alignment[] donorAlign, boolean print) {
        OpenBitSet[][] result=new OpenBitSet[donorAlign.length][4];
        for (int da = 0; da < result.length; da++) {
            int donorOffset=unimpAlign.getSiteOfPhysicalPosition(donorAlign[da].getPositionInLocus(0), donorAlign[da].getLocus(0), donorAlign[da].getSNPID(0));
            OpenBitSet goodMask=new OpenBitSet(donorAlign[da].getSiteCount());
            OpenBitSet swapMjMnMask=new OpenBitSet(donorAlign[da].getSiteCount());
            OpenBitSet errorMask=new OpenBitSet(donorAlign[da].getSiteCount());
            OpenBitSet invariantMask=new OpenBitSet(donorAlign[da].getSiteCount());
            int siteConflicts=0, swaps=0, invariant=0, good=0;
            for (int i = 0; i < donorAlign[da].getSiteCount(); i++) {
                /*we have three classes of data:  invariant in one alignment, conflicts about minor and minor,
                *swaps of major and minor.  Adding the invariant reduces imputation accuracy.
                *the major/minor swaps should be flipped in the comparisons
                */
                byte tMj=unimpAlign.getMajorAllele(i+donorOffset);
                byte tMn=unimpAlign.getMinorAllele(i+donorOffset);
                byte daMj=donorAlign[da].getMajorAllele(i);
                byte daMn=donorAlign[da].getMinorAllele(i);
                if(daMn==Alignment.UNKNOWN_ALLELE) {
                    invariant++;
                    invariantMask.set(i);
                    goodMask.set(i);
                } else 
                if((daMj==tMn)&&(daMn==tMj)) {
                    swaps++;
                    swapMjMnMask.set(i);
                    goodMask.set(i);
                } else
                if((daMj!=tMj)) {
                    siteConflicts++;
                    errorMask.set(i);
                    goodMask.set(i);
                } 

            }
            goodMask.not();
            if(print) System.out.println("Donor:"+da+"invariant in donor:"+invariant+" swapConflicts:"+swaps+" errors:"+siteConflicts);
            result[da]=new OpenBitSet[] {goodMask, swapMjMnMask, errorMask, invariantMask};
        }
        return result;
    }
    
    
    
    private int[] countUnknownAndHets(byte[] a) {
        int cnt=0, cntHets=0;
        for (int i = 0; i < a.length; i++) {
            if(a[i]==Alignment.UNKNOWN_DIPLOID_ALLELE) {cnt++;}
            else if(AlignmentUtils.isHeterozygous(a[i])) {cntHets++;}
        }
        return new int[]{cnt,cntHets};
    }
    
    private BitSet[] arrangeMajorMinorBtwAlignments(Alignment unimpAlign, int bt, int donorOffset, int donorLength,
            OpenBitSet goodMask, OpenBitSet swapMjMnMask) { 
        int unimpAlignStartBlock=donorOffset/64;
        int shift=(donorOffset-(unimpAlignStartBlock*64));
        int unimpAlignEndBlock=unimpAlignStartBlock+((donorLength+shift-1)/64);
        OpenBitSet mjUnImp=new OpenBitSet(unimpAlign.getAllelePresenceForSitesBlock(bt, 0, unimpAlignStartBlock, unimpAlignEndBlock+1));
        OpenBitSet mnUnImp=new OpenBitSet(unimpAlign.getAllelePresenceForSitesBlock(bt, 1, unimpAlignStartBlock, unimpAlignEndBlock+1));
        OpenBitSet mjTbs=new OpenBitSet(donorLength);
        OpenBitSet mnTbs=new OpenBitSet(donorLength);
        for (int i = 0; i < donorLength; i++) {
            if(mjUnImp.fastGet(i+shift)) mjTbs.set(i);
            if(mnUnImp.fastGet(i+shift)) mnTbs.set(i);
        }
        OpenBitSet newmj=new OpenBitSet(mjTbs);
        OpenBitSet newmn=new OpenBitSet(mnTbs);
        mjTbs.and(goodMask);
        mnTbs.and(goodMask);
 //       System.out.printf("mjTbs:%d Goodmask:%d swapMjMnMask:%d",mjTbs.getNumWords(),goodMask.getNumWords(), swapMjMnMask.getNumWords());
        if(isSwapMajorMinor) {
            newmj.and(swapMjMnMask);
            newmn.and(swapMjMnMask);
            mjTbs.or(newmn);
            mnTbs.or(newmj);
        }
//        System.out.printf("Arrange shift:%d finalEnd:%d totalBits:%d DesiredLength:%d %n", shift, finalEnd, totalBits, donorLength);
        BitSet[] result={mjTbs,mnTbs};
        return result;
    }
 
    /**
     * If the target regions has Mendelian errors that it looks for overlapping regional
     * solutions that are better.
     * TODO:  This can be done much more robustly.
     * @param impT
     * @param targetTaxon
     * @param regionHypth 
     */ 
    private ImputedTaxon applyBlockNN(ImputedTaxon impT, int targetTaxon, 
            Alignment donorAlign, int donorOffset, DonorHypoth[][] regionHypth, boolean hybridMode, BitSet[] maskedTargetBits, 
            int minMinorCnt, double maxHybridErrorRate, int[] donorIndices) {
        int focusBlockdone=0, inbred=0, hybrid=0, noData=0;
        int blocks=maskedTargetBits[0].getNumWords();
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            if((regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()<maximumInbredError)) {
                setAlignmentWithDonors(donorAlign, regionHypth[focusBlock],donorOffset, true,impT);
                focusBlockdone++;
                impT.incBlocksSolved();
                inbred++;
            }  else {
                int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),
                        maskedTargetBits[1].getBits(), focusBlock, minMinorCnt);
                if(resultRange==null) {noData++; continue; }//no data in the focus Block
                 //search for the best hybrid donors for a segment
                if(hybridMode&&(regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()>maximumInbredError)) {
                    long[] mjTRange=maskedTargetBits[0].getBits(resultRange[0],resultRange[2]);
                    long[] mnTRange=maskedTargetBits[1].getBits(resultRange[0],resultRange[2]);
                    regionHypth[focusBlock]=getBestHybridDonors(targetTaxon, mjTRange, mnTRange, resultRange[0],resultRange[2], 
                            focusBlock, donorAlign, regionHypth[focusBlock], false, donorIndices);
                }
          //      if((regionHypth[focusBlock][0]!=null)) System.out.printf("%d %d %g %n",targetTaxon, focusBlock, regionHypth[focusBlock][0].getErrorRate() );
                if((regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()<maxHybridErrorRate)) {
                    setAlignmentWithDonors(donorAlign, regionHypth[focusBlock], donorOffset, true,impT);
                    focusBlockdone++;
                    impT.incBlocksSolved();
                    hybrid++;
                } else {
                    if(regionHypth[focusBlock][0]==null) {noData++;} 
//                    else{
//                        System.out.println(regionHypth[focusBlock][0].getErrorRate());
//                    }
                }
            }
        }
        
        int leftNullCnt=0;
        for (DonorHypoth[] dh : regionHypth) {if(dh[0]==null) leftNullCnt++;}
        if(testing==1) System.out.printf("targetTaxon:%d hybridError:%g block:%d focusBlockdone:%d null:%d inbredDone:%d hybridDone:%d noData:%d %n",
                targetTaxon, maxHybridErrorRate, blocks,focusBlockdone, leftNullCnt, inbred, hybrid, noData);
        return impT;
    }
    
    private DonorHypoth getStateBasedOnViterbi(DonorHypoth dh, int donorOffset, Alignment donorAlign) {
        TransitionProbability tp = new TransitionProbability();
        EmissionProbability ep = new EmissionProbability();
        tp.setTransitionProbability(transition);
        ep.setEmissionProbability(emission);
        int startSite=dh.startBlock*64;
        int endSite=(dh.endBlock*64)+63;
        if(endSite>=donorAlign.getSiteCount()) endSite=donorAlign.getSiteCount()-1;
        int sites=endSite-startSite+1;
        byte[] calls=new byte[sites];
        //System.out.printf("%d %d %d %n",dh.donor1Taxon, startSite, endSite+1);
        byte[] d1b=donorAlign.getBaseRange(dh.donor1Taxon, startSite, endSite+1);
        byte[] d2b=donorAlign.getBaseRange(dh.donor2Taxon, startSite, endSite+1);
        byte[] t1b=unimpAlign.getBaseRange(dh.targetTaxon, startSite+donorOffset, endSite+1+donorOffset);
        int informSites=0, nonMendel=0;
        ArrayList<Byte> nonMissingObs = new ArrayList<Byte>();
        ArrayList<Integer> snpPositions = new ArrayList<Integer>();
        for(int cs=0; cs<sites; cs++) {
            if(t1b[cs]==Alignment.UNKNOWN_DIPLOID_ALLELE) continue; 
            if(d1b[cs]==Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
            if(d2b[cs]==Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
            if(t1b[cs]==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) continue; 
            if(d1b[cs]==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) continue;
            if(d2b[cs]==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) continue;
            if(d1b[cs]==d2b[cs]) {
                if(t1b[cs]!=d1b[cs]) nonMendel++;
                continue;
            }
            informSites++;
            byte state=1;
            if(t1b[cs]==d1b[cs]) {state=0;}
            else if(t1b[cs]==d2b[cs]) {state=2;}
            nonMissingObs.add(state);
            snpPositions.add(cs+startSite);
        }
        if(informSites<10) return null;
        double nonMendRate=(double)nonMendel/(double)informSites;
        if(testing==1) System.out.printf("NonMendel:%d InformSites:%d ErrorRate:%g %n",nonMendel, informSites, nonMendRate);
        if(nonMendRate>this.maximumInbredError*5) return null;
        byte[] informStates=new byte[informSites];
        for (int i = 0; i < informStates.length; i++) informStates[i]=nonMissingObs.get(i);
        int[] pos=new int[informSites];
        for (int i = 0; i < pos.length; i++) pos[i]=snpPositions.get(i);
        int chrlength = donorAlign.getPositionInLocus(endSite) - donorAlign.getPositionInLocus(startSite);
        tp.setAverageSegmentLength( chrlength / sites );
        tp.setPositions(pos);
        
	double probHeterozygous=0.5;
        double phom = (1 - probHeterozygous) / 2;
        double[] pTrue = new double[]{phom, .25*probHeterozygous ,.5 * probHeterozygous, .25*probHeterozygous, phom};
		
        ViterbiAlgorithm va = new ViterbiAlgorithm(informStates, tp, ep, pTrue);
	va.calculate();
        if(testing==1) System.out.println("Input:"+Arrays.toString(informStates));
        byte[] resultStates=va.getMostProbableStateSequence();
        if(testing==1) System.out.println("Resul:"+Arrays.toString(resultStates));
        DonorHypoth dh2=new DonorHypoth(dh.targetTaxon,dh.donor1Taxon, 
                        dh.donor2Taxon, dh.startBlock, dh.focusBlock, dh.endBlock);
        int currPos=0;
        for(int cs=0; cs<sites; cs++) {
            calls[cs]=(resultStates[currPos]==1)?(byte)1:(byte)(resultStates[currPos]/2); //converts the scale back to 0,1,2 from 0..4
            if((pos[currPos]<cs+startSite)&&(currPos<resultStates.length-1)) currPos++;
        }
        dh2.phasedResults=calls;
        return dh2;
    }
    
    
    
    /**
     * Given a start 64 site block, it expands to the left and right until it hits
     * the minimum Minor Site count in the target taxon
     * @param mnT - minor allele bit presence in a series of longs
     * @param focusBlock
     * @param minMinorCnt
     * @return arrays of blocks {startBlock, focusBlock, endBlock}
     */
    private int[] getBlockWithMinMinorCount(long[] mjT, long[] mnT, int focusBlock, int minMinorCnt) {
        int blocks=mjT.length;
        int majorCnt=Long.bitCount(mjT[focusBlock]);
        int minorCnt=Long.bitCount(mnT[focusBlock]);
        int endBlock=focusBlock, startBlock=focusBlock;
        int minMajorCnt=minMinorCnt*minMajorRatioToMinorCnt;
        while((minorCnt<minMinorCnt)&&(majorCnt<minMajorCnt)) {
            boolean preferMoveStart=(focusBlock-startBlock<endBlock-focusBlock)?true:false;
            if(startBlock==0) {preferMoveStart=false;}
            if(endBlock==blocks-1) {preferMoveStart=true;}
            if((startBlock==0)&&(endBlock==blocks-1)) break;
            if(preferMoveStart) {//expand start
                startBlock--;
                minorCnt+=Long.bitCount(mnT[startBlock]);
                majorCnt+=Long.bitCount(mjT[startBlock]);
            } else { //expand end
                endBlock++;
                minorCnt+=Long.bitCount(mnT[endBlock]); 
                majorCnt+=Long.bitCount(mjT[startBlock]);
            } 
        }
        int[] result={startBlock, focusBlock, endBlock};
        return result;
    }
  
    
    /**
     * Calculates an inbred distance
     * TODO:  Should be converted to the regular same, diff, hets, & sum
     * @param modBitsOfTarget
     * @return array with sites, same count, diff count, het count
     */
    private void calcInbredDist(ImputedTaxon impT, BitSet[] modBitsOfTarget, Alignment donorAlign) {
        int blocks=modBitsOfTarget[0].getNumWords();
        impT.allDist=new byte[donorAlign.getSequenceCount()][4][blocks];
        long[] iMj=modBitsOfTarget[0].getBits();
        long[] iMn=modBitsOfTarget[1].getBits();
        for (int donor1 = 0; donor1 < impT.allDist.length; donor1++) {
            long[] jMj=donorAlign.getAllelePresenceForAllSites(donor1, 0).getBits();
            long[] jMn=donorAlign.getAllelePresenceForAllSites(donor1, 1).getBits();
            for (int i = 0; i <blocks; i++) {
                long same = (iMj[i] & jMj[i]) | (iMn[i] & jMn[i]);
                long diff = (iMj[i] & jMn[i]) | (iMn[i] & jMj[i]);
                long hets = same & diff;
                int sameCnt = BitUtil.pop(same);
                int diffCnt = BitUtil.pop(diff);
                int hetCnt = BitUtil.pop(hets);
                int sites = sameCnt + diffCnt - hetCnt;
                impT.allDist[donor1][2][i]=(byte)diffCnt;
                impT.allDist[donor1][3][i]=(byte)hetCnt;
                impT.allDist[donor1][0][i]=(byte)sites;
                impT.allDist[donor1][1][i]=(byte)sameCnt;
            }
        }
    }
    
    /**
     * Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}
     */
    private DonorHypoth[] getBestInbredDonors(int targetTaxon, ImputedTaxon impT, int startBlock, int endBlock, 
            int focusBlock, Alignment donorAlign, int[] donor1indices) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        int donorTaxaCnt=donorAlign.getSequenceCount();
        for (int d1 : donor1indices) {
            int testSites=0;
            int sameCnt = 0, diffCnt = 0, hetCnt = 0;
            for (int i = startBlock; i <=endBlock; i++) {
                sameCnt+=impT.allDist[d1][1][i];
                diffCnt+=impT.allDist[d1][2][i];
                hetCnt+=impT.allDist[d1][3][i];
            }
            testSites= sameCnt + diffCnt - hetCnt;
            if(testSites<minTestSites) continue;     
            double testPropUnmatched = 1.0-(((double) (sameCnt) - (double)(0.5*hetCnt)) / (double) (testSites));
            int totalMendelianErrors=(int)((double)testSites*testPropUnmatched);
            inc+=1e-9;  //this looks strange, but makes the keys unique and ordered
            testPropUnmatched+=inc;  //this looks strange, but makes the keys unique and ordered 
            if(testPropUnmatched<lastKeytestPropUnmatched) {
                DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, startBlock,
                    focusBlock, endBlock, testSites, totalMendelianErrors);
                DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                if(prev!=null) System.out.println("Inbred TreeMap index crash:"+testPropUnmatched);
                if(bestDonors.size()>maxDonorHypotheses) {
                    bestDonors.remove(bestDonors.lastKey());
                    lastKeytestPropUnmatched=bestDonors.lastKey();
                }
            }  
        }
        //TODO consider reranking by Poisson.  The calculating Poisson on the fly is too slow.
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
        return result;
    }
    
    /**
     * Produces a sort list of most prevalent donorHypotheses across the donorAlign.
     * Currently it is only looking at the very best for each focus block.
     * @param allDH
     * @param minHypotheses
     * @return 
     */
    private int[] getAllBestDonorsAcrossChromosome(DonorHypoth[][] allDH, int minHypotheses) {
        TreeMap<Integer,Integer> bd=new TreeMap<Integer,Integer>();
        for (int i = 0; i < allDH.length; i++) {
            DonorHypoth dh=allDH[i][0];
            if(dh==null) continue;
            if(bd.containsKey(dh.donor1Taxon)) {
                bd.put(dh.donor1Taxon, bd.get(dh.donor1Taxon)+1);
            } else {
                bd.put(dh.donor1Taxon, 1);
            }
        }
        if(testing==1) System.out.println(bd.size()+":"+bd.toString());
        if(testing==1) System.out.println("");
        int highDonors=0;
        for (int i : bd.values()) {if(i>minHypotheses) highDonors++;}
        int[] result=new int[highDonors];    
        highDonors=0;
        for (Map.Entry<Integer,Integer> e: bd.entrySet()) {
            if(e.getValue()>minHypotheses) result[highDonors++]=e.getKey();}
        if(testing==1) System.out.println(Arrays.toString(result));
        return result;
    }
  
    private DonorHypoth[] getBestHybridDonors(int targetTaxon, long[] mjT, long[] mnT, 
            int startBlock, int endBlock, int focusBlock, Alignment donorAlign, DonorHypoth[] inbredHypoth, 
            boolean compareDonorDist, int[] donorIndices) {
        int nonNull=0;
        for (DonorHypoth dH : inbredHypoth) {
            if(dH!=null) nonNull++;
        }
        int[] inbredIndices=new int[inbredHypoth.length];
        for (int i = 0; i < nonNull; i++) {
            inbredIndices[i]=inbredHypoth[i].donor1Taxon;   
        }
        return getBestHybridDonors(targetTaxon, mjT, mnT, startBlock,  endBlock, focusBlock, donorAlign,
                inbredIndices, donorIndices, compareDonorDist);
    }
    
    /**
     * Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}
     */
    private DonorHypoth[] getBestHybridDonors(int targetTaxon, long[] mjT, long[] mnT, 
            int startBlock, int endBlock, int focusBlock, Alignment donorAlign, int[] donor1Indices, int[] donor2Indices,
            boolean viterbiSearch) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        double[] donorDist;
        for (int d1: donor1Indices) {
            long[] mj1=donorAlign.getAllelePresenceForSitesBlock(d1, 0, startBlock, endBlock+1);
            long[] mn1=donorAlign.getAllelePresenceForSitesBlock(d1, 1, startBlock, endBlock+1);
            for (int d2 : donor2Indices) {
                if((!viterbiSearch)&&(d1==d2)) continue;
                long[] mj2=donorAlign.getAllelePresenceForSitesBlock(d2, 0, startBlock, endBlock+1);
                long[] mn2=donorAlign.getAllelePresenceForSitesBlock(d2, 1, startBlock, endBlock+1);
                if(viterbiSearch) {
                    donorDist=IBSDistanceMatrix.computeHetBitDistances(mj1, mn1, mj2, mn2, minTestSites);
                    if((d1!=d2)&&(donorDist[0]<this.maximumInbredError)) continue;
                } 
                int[] mendErr=mendelErrorComparison(mjT, mnT, mj1, mn1, mj2, mn2);        
                if(mendErr[1]<minTestSites) continue;
                double testPropUnmatched=(double)(mendErr[0])/(double)mendErr[1];
                inc+=1e-9;
                testPropUnmatched+=inc; 
                if(testPropUnmatched<lastKeytestPropUnmatched) {
                    DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, startBlock,
                        focusBlock, endBlock, mendErr[1], mendErr[0]);
                    DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                    if(prev!=null) {
                        System.out.println("Hybrid TreeMap index crash:"+testPropUnmatched);
                    }
                    if(bestDonors.size()>maxDonorHypotheses) {
                        bestDonors.remove(bestDonors.lastKey());
                        lastKeytestPropUnmatched=bestDonors.lastKey();
                    }
                }  
            }
        }
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
        return result;
    }
    
    private int[] mendelErrorComparison(long[] mjT, long[] mnT, long[] mj1, long[] mn1, 
            long[] mj2, long[] mn2) {
        int mjUnmatched=0;
        int mnUnmatched=0;
        int testSites=0;
        for (int i = 0; i < mjT.length; i++) {
            long siteMask=(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
            mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
            mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
            testSites+=Long.bitCount(siteMask);
        }
        int totalMendelianErrors=mjUnmatched+mnUnmatched;
       // double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
        return (new int[] {totalMendelianErrors, testSites});       
    }
    
    
    private ImputedTaxon setAlignmentWithDonors(Alignment donorAlign, DonorHypoth[] theDH, int donorOffset,
            boolean setJustFocus, ImputedTaxon impT) {
        if(theDH[0].targetTaxon<0) return impT;
        boolean print=false;
        int startSite=(setJustFocus)?theDH[0].getFocusStartSite():theDH[0].startSite;
        int endSite=(setJustFocus)?theDH[0].getFocusEndSite():theDH[0].endSite;
        if(endSite>=donorAlign.getSiteCount()) endSite=donorAlign.getSiteCount()-1;
//        if (print) System.out.println("B:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
        int[] prevDonors;
        if(impT.breakPoints.isEmpty()) {
            prevDonors=new int[]{-1, -1};
        } else {
            prevDonors=impT.breakPoints.lastEntry().getValue();
        }
        int[] currDonors=prevDonors;       
        for(int cs=startSite; cs<=endSite; cs++) {
            byte donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE; 
            byte neighbor=0;
            for (int i = 0; (i < theDH.length) && (donorEst==Alignment.UNKNOWN_DIPLOID_ALLELE); i++) {
                neighbor++;
                if((theDH[i]==null)||(theDH[i].donor1Taxon<0)) continue;
                if(theDH[i].getErrorRate()>this.maximumInbredError) continue;
                byte bD1=donorAlign.getBase(theDH[i].donor1Taxon, cs);
                if(theDH[i].getPhaseForSite(cs)==0) {
                    donorEst=bD1;
                    if(i==0) currDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor1Taxon};
                }
                else {
                    byte bD2=donorAlign.getBase(theDH[i].donor2Taxon, cs);
                    if(theDH[i].getPhaseForSite(cs)==2) {
                        donorEst=bD2;
                        if(i==0) currDonors=new int[]{theDH[0].donor2Taxon, theDH[0].donor2Taxon};
                    }
                    else {
                        donorEst=AlignmentUtils.getUnphasedDiploidValueNoHets(bD1, bD2);
                        if(i==0) currDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor2Taxon};
                    }
                }
            }
            byte knownBase=impT.getOrigGeno(cs+donorOffset);
            if(!Arrays.equals(prevDonors, currDonors)) {
                impT.breakPoints.put(donorAlign.getPositionInLocus(cs), currDonors);
                prevDonors=currDonors;
            }
            if(theDH[0].phasedResults==null) {impT.chgHis[cs+donorOffset]=(byte)-neighbor;}
            else {impT.chgHis[cs+donorOffset]=(byte)neighbor;}
            
            impT.impGeno[cs+donorOffset]= donorEst;  //predicted based on neighbor
            if(knownBase==Alignment.UNKNOWN_DIPLOID_ALLELE) {impT.resolveGeno[cs+donorOffset]= donorEst;}
            else {if(AlignmentUtils.isHeterozygous(donorEst)) {
                if(resolveHetIfUndercalled&&AlignmentUtils.isPartiallyEqual(knownBase,donorEst)) 
                {//System.out.println(theDH[0].targetTaxon+":"+knownBase+":"+donorEst);
                    impT.resolveGeno[cs+donorOffset]= donorEst;}
            }}
        } //end of cs loop
        //enter a stop of the DH at the beginning of the next block
        int lastDApos=donorAlign.getPositionInLocus(endSite);
        int nextSite=unimpAlign.getSiteOfPhysicalPosition(lastDApos, donorAlign.getLocus(0))+1;
        if(nextSite<unimpAlign.getSiteCount()) impT.breakPoints.put(unimpAlign.getPositionInLocus(nextSite), new int[]{-1,-1});
    //    if (print) System.out.println("E:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
        return impT;
    }
    
    private double calcErrorForTaxonAndSite(ImputedTaxon impT) {
        for (int cs = 0; cs < impT.getOrigGeno().length; cs++) {
            byte donorEst=impT.getImpGeno(cs);
            byte knownBase=impT.getOrigGeno(cs);
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                if(AlignmentUtils.isHeterozygous(donorEst)||AlignmentUtils.isHeterozygous(knownBase)) {
                    totalHets++;
                } else if(knownBase==donorEst) {
                    totalRight++;
                    siteCorrectCnt[cs]++;
                    taxonCorrectCnt[impT.taxon]++;
                } else {
                    totalWrong++;
                    siteErrors[cs]++;
                    taxonErrors[impT.taxon]++;
                }
            }
            
        }
        return (double)taxonErrors[impT.taxon]/(double)(taxonCorrectCnt[impT.taxon]+taxonErrors[impT.taxon]);
    }
    
    public static int[] compareAlignment(String origFile, String maskFile, String impFile, boolean noMask) {
        boolean taxaOut=false;
        Alignment oA=ImportUtils.readGuessFormat(origFile, false);
        System.out.printf("Orig taxa:%d sites:%d %n",oA.getSequenceCount(),oA.getSiteCount());        
        Alignment mA=null;
        if(noMask==false) {mA=ImportUtils.readGuessFormat(maskFile, false);
            System.out.printf("Mask taxa:%d sites:%d %n",mA.getSequenceCount(),mA.getSiteCount());
        }
        Alignment iA=ImportUtils.readGuessFormat(impFile, false);
        System.out.printf("Imp taxa:%d sites:%d %n",iA.getSequenceCount(),iA.getSiteCount());
        int correct=0;
        int errors=0;
        int unimp=0;
        int hets=0;
        int gaps=0;
        for (int t = 0; t < iA.getSequenceCount(); t++) {
            int e=0,c=0,u=0,h=0;
            int oATaxa=oA.getIdGroup().whichIdNumber(iA.getFullTaxaName(t));
            for (int s = 0; s < iA.getSiteCount(); s++) {
                if(noMask||(oA.getBase(oATaxa, s)!=mA.getBase(t, s))) {
                    byte ib=iA.getBase(t, s);
                    byte ob=oA.getBase(oATaxa, s);
                    if((ib==Alignment.UNKNOWN_DIPLOID_ALLELE)||(ob==Alignment.UNKNOWN_DIPLOID_ALLELE)) {unimp++; u++;}
                    else if(ib==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {gaps++;}
                    else if(ib==ob) {
                        correct++;
                        c++;
                    } else {
                        if(AlignmentUtils.isHeterozygous(ob)||AlignmentUtils.isHeterozygous(ib)) {hets++; h++;}
                        else {errors++; 
                            e++;
//                            if(t==0) System.out.printf("%d %d %s %s %n",t,s,oA.getBaseAsString(oATaxa, s), iA.getBaseAsString(t, s));
                        }
                    }
                }       
            }
            if(taxaOut) System.out.printf("%s %d %d %d %d %n",iA.getTaxaName(t),u,h,c,e);
        }
        System.out.println("MFile\tIFile\tGap\tUnimp\tUnimpHets\tCorrect\tErrors");   
        System.out.printf("%s\t%s\t%d\t%d\t%d\t%d\t%d%n",maskFile, impFile, gaps, unimp,hets,correct,errors);   
        return new int[]{gaps, unimp,hets,correct,errors};
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-d", "--donorH", true);
        engine.add("-sC", "--startChrom", true);
        engine.add("-eC", "--endChrom", true);
        engine.add("-minMnCnt", "--minMnCnt", true);
        engine.add("-mxInbErr", "--mxInbErr", true);
        engine.add("-mxHybErr", "--mxHybErr", true);
        engine.add("-inbNNOff", "--inbNNOff", false);
        engine.add("-hybNNOff", "--hybNNOff", false);
        engine.add("-mxDonH", "--mxDonH", true);
        engine.add("-mnTestSite", "--mnTestSite", true);
        engine.add("-projA", "--projAlign", false);
        engine.add("-runChrMode", "--runChrMode", false);
        engine.parse(args);
        hmpFile = engine.getString("-hmp");
        outFileBase = engine.getString("-o");
        donorFile = engine.getString("-d");
        if (engine.getBoolean("-sC")) {
            startChr = Integer.parseInt(engine.getString("-sC"));
        }
        if (engine.getBoolean("-eC")) {
            endChr = Integer.parseInt(engine.getString("-eC"));
        }
        if (engine.getBoolean("-mxInbErr")) {
            maximumInbredError = Double.parseDouble(engine.getString("-mxInbErr"));
        }
        if (engine.getBoolean("-mxHybErr")) {
            maxHybridErrorRate = Double.parseDouble(engine.getString("-mxHybErr"));
        }
        if (engine.getBoolean("-minMnCnt")) {
            minMinorCnt = Integer.parseInt(engine.getString("-minMnCnt"));
        }
        if (engine.getBoolean("-inbNNOff")) inbredNN=false;
        if (engine.getBoolean("-hybNNOff")) hybridNN=false;
        if (engine.getBoolean("-mxDonH")) {
            maxDonorHypotheses = Integer.parseInt(engine.getString("-mxDonH"));
        }
        if (engine.getBoolean("-mnTestSite")) {
            minTestSites = Integer.parseInt(engine.getString("-mnTestSite"));
        }
        if (engine.getBoolean("-projA")) isOutputProjection=true;
        if (engine.getBoolean("-runChrMode")) isChromosomeMode=true;
    }



    private void printUsage() {
        myLogger.info(
                "\n\n\nAvailable options for the BiParentalErrorCorrectionPlugin are as follows:\n"
                + "-hmp   Input HapMap file(s) 'c+' or 'gX' to denote variable chromosomes\n"
                + "-d    Donor haplotype files 'c+s+' or 'gX' to denote sections\n"
                + "-o     Output HapMap file(s) 'c+' to denote variable chromosomes\n"
                + "-runChrMode Deprecated mode to run individual chromosomes\n"
                + "-sC    Start chromosome\n"
                + "-eC    End chromosome\n"
                + "-minMnCnt    Minor number of minor alleles in the search window (or "+minMajorRatioToMinorCnt+"X major)\n"
                + "-mxInbErr    Maximum inbred error rate\n"
                + "-mxHybErr    Maximum hybrid error rate\n"
                + "-inbNNOff    Whether to use inbred NN (default:"+inbredNN+")\n"
                + "-hybNNOff    Whether to use both the hybrid NN (default:"+hybridNN+")\n"
                + "-mxDonH   Maximum number of donor hypotheses to be explored (default: "+maxDonorHypotheses+")\n"
                + "-mnTestSite   Minimum number of sites to test for NN IBD (default:"+minTestSites+")\n"
                + "-projA   Create a projection alignment for high density markers (default off)\n"
                );
    }

    @Override
    public DataSet performFunction(DataSet input) {
        if(isChromosomeMode) {
            for (int chr = startChr; chr <=endChr; chr++) {
               String chrHmpFile=hmpFile.replace("chr+", "chr"+chr);
               chrHmpFile=chrHmpFile.replace("c+", "c"+chr);
               String chrDonorFile=donorFile.replace("chr+", "chr"+chr);
               chrDonorFile=chrDonorFile.replace("c+", "c"+chr);
               String chrOutFile=outFileBase.replace("chr+", "chr"+chr);
               chrOutFile=chrOutFile.replace("c+", "c"+chr);

               runMinorWindowViterbiImputation(chrDonorFile, chrHmpFile, chrOutFile, 
                       minMinorCnt, minTestSites, 100, maxHybridErrorRate, isOutputProjection, false);
            }
        } else {
            runMinorWindowViterbiImputation(donorFile, hmpFile, outFileBase, 
                       minMinorCnt, minTestSites, 100, maxHybridErrorRate, isOutputProjection, false);
        }
        return null;
    }
    


    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "ImputeByNN&HMM";
    }

    @Override
    public String getToolTipText() {
        return "Imputation that relies on a combination of HMM and Nearest Neighbor";
    }
    
    public static void main(String[] args) {
        
//        String rootOrig="/Volumes/LaCie/build20120701/IMP26/orig/";
//        String rootHaplos="/Volumes/LaCie/build20120701/IMP26/haplos/";
//        String rootImp="/Volumes/LaCie/build20120701/IMP26/imp/";
        String rootOrig="/Users/edbuckler/SolexaAnal/GBS/build20120701/IMP26/orig/";
        String rootHaplos="/Users/edbuckler/SolexaAnal/GBS/build20120701/IMP26/haplos/";
        String rootImp="/Users/edbuckler/SolexaAnal/GBS/build20120701/IMP26/imp/";
        //String unImpTargetFile=rootOrig+"AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr+.hmp.txt.gz";
  //      String unImpTargetFile=rootOrig+"Samp82v26.chr8.hmp.txt.gz";
//        String unImpTargetFile=rootOrig+"AllZeaGBS_v2.6.chr+.hmp.h5";
       String unImpTargetFile=rootOrig+"USNAM142v26.chr10.hmp.txt.gz";
//       String unImpTargetFile=rootOrig+"Ames105v26.chr10.hmp.txt.gz";
        String donorFile=rootHaplos+"all26_8k.c10s+.hmp.txt.gz";
//        String donorFile=rootHaplos+"HM26_Allk.c10s+.hmp.txt.gz";
//        String impTargetFile=rootImp+"newmnaTall26.c+.imp.hmp.txt.gz";
       Random r=new Random();
       String impTargetFile=rootImp+r.nextInt()+"newmnaTall26.c+.imp.hmp.h5";
  //      String impTargetFile=rootImp+"T3AllZeaGBSv2_6.c+.pa.txt.gz";
        
      
        String[] args2 = new String[]{
            "-hmp", unImpTargetFile,
            "-d",donorFile,
            "-o", impTargetFile,
            "-sC","10",
            "-eC","10",
            "-minMnCnt","20",
            "-mxInbErr","0.02",
            "-mxHybErr","0.005",
            "-mnTestSite","50",
            "-mxDonH","10",
            "-runChrMode",
 //           "-projA",
        };
        System.out.println("TasselPrefs:"+TasselPrefs.getAlignmentRetainRareAlleles());
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        MinorWindowViterbiImputationPlugin plugin = new MinorWindowViterbiImputationPlugin();
        plugin.setParameters(args2);
        plugin.performFunction(null);
 //       compareAlignment(rootImp+"newmnaTall26.c10.imp.hmp.txt.gz", null, rootImp+"newmnaTall26.c10.imp.mhmp.h5", true);
//        TasselPrefs.putAlignmentRetainRareAlleles(false);
//        System.out.println("Reading PA file");
////        ProjectionAlignment pa=ProjectionAlignment.getInstance(rootImp+"TAllZeaGBSv2_6.c10.pa.txt.gz", donorFile);
//        ProjectionAlignment pa=ProjectionAlignment.getInstance(rootImp+"AllZeaGBSv2_6.c10.pa.txt.gz", donorFile);
//        System.out.println("Writing PA file to HapMap");
//        ExportUtils.writeToHapmap(pa, false, rootImp+"PostPATall26.c10.hmp.txt.gz", '\t', null);
//        System.out.printf("Cache:%d Lookup:%d %n",pa.cacheUseCnt, pa.lookupCnt);
//        compareAlignment(unImpTargetFile, null, rootImp+"PostPATall26.c10.hmp.txt.gz", true);
    }
}

class ImputedTaxon {
    int taxon;
    private boolean isProjection;
    private boolean segmentSolved=false; //change to true when done
    private byte[] origGeno;
    byte[] impGeno;  //the best imputed estimate, orginal sequence is not part of this
    byte[] resolveGeno;  //what to set the alignment with combination of original and imp.
    byte[] chgHis; //Viterbi is negative, blockNN is positive
    BitSet[] modBitsOfTarget;
    byte[][][] allDist;
     TreeMap<Integer,int[]> breakPoints;
     private int blocksSolved=0; 

    public ImputedTaxon(int taxon, byte[] origGeno, boolean isProjection) {
        this.taxon = taxon;
        this.origGeno = origGeno;
        this.isProjection=isProjection;
        impGeno=Arrays.copyOf(origGeno, origGeno.length); //imputed sequence
        resolveGeno=Arrays.copyOf(origGeno, origGeno.length); //imputed sequence
        chgHis=new byte[origGeno.length];
        breakPoints=new TreeMap<Integer,int[]>();
        breakPoints.put(0, new int[]{-1,-1});
    }

    public byte[] getOrigGeno() {
        return origGeno;
    }
    
    public byte getOrigGeno(int site) {
        return origGeno[site];
    }

    public byte[] getImpGeno() {
        return impGeno;
    }
    
    public byte getImpGeno(int site) {
        return impGeno[site];
    }

    public void setImpGeno(int site, byte diploidGenotype) {
        impGeno[site] = diploidGenotype;
    }
    
    public boolean isProjection() {
        return isProjection;
    }

    public boolean isSegmentSolved() {
        return segmentSolved;
    }

    public void setSegmentSolved(boolean segmentSolved) {
        this.segmentSolved = segmentSolved;
    }

    public int getBlocksSolved() {
        return blocksSolved;
    }
    
    public void incBlocksSolved() {
        blocksSolved++;
    }
    
    
    
}
