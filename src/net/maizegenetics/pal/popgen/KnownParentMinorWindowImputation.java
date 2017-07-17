package net.maizegenetics.pal.popgen;

import java.util.ArrayList;
import java.util.Random;
import java.util.TreeMap;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;


/**
 *  Imputation methods that relies on a list of possible parents, and optimized for
 * highly heterozygous samples.  It works at the scale of 64 sites to accelerate searching
 * through all possible parental combinations for a window.  
 * 
 * It only begins each 64 site block and expand outward to find set number of minor
 * alleles in the target sequence, and then it looks for all possible parents.
 * 
 * @deprecated Use MinorWindowViterbiImputationPlugin {@link net.maizegenetics.gbs.pipeline.MinorWindowViterbiImputationPlugin}
 * @author edbuckler, kswarts, aromero
 */
@Deprecated
public class KnownParentMinorWindowImputation {
    private Alignment unimpAlign;
    private Alignment donorAlign;
    private int testing=0;  //level of reporting to stdout
    private OpenBitSet swapMjMnMask;  //mask of sites that major and minor are swapped
    private int parentsRight=0, parentsWrong=0;
    private int blocks=-1;
    int[] hSite, hTaxon;
    byte[] hState;
    int maskSitCnt=0;
    boolean maskAndTest=true;
    private double maximumInbredError=0.02;  //inbreds are tested first, if too much error hybrids are tested.
    
    private static int highMask=0xF0;
    private static int lowMask=0x0F;

    /**
     * 
     * @param donorFile should be imputed inbreds that were founders of population
     * @param unImpTargetFile sites must match exactly with donor file
     * @param exportFile 
     * @param minMinorCnt determines the size of the search window, low recombination 20-30, high recombination 10-15
     * @param resolveMethod 0=focus sites; 1=regional solution (better)
     */
    public KnownParentMinorWindowImputation(String donorFile, String unImpTargetFile, 
            String exportFile, int minMinorCnt, int resolveMethod) {
        donorAlign=ImportUtils.readFromHapmap(donorFile, false, (ProgressListener)null);
        donorAlign.optimizeForTaxa(null);
        unimpAlign=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
        unimpAlign.optimizeForTaxa(null);
        
        if(maskAndTest) maskSites(300);  //mask sites to test imputation accuracy
        blocks=unimpAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
        
        //Create mask for all sites where major & minor are swapped in alignments
        swapMjMnMask=new OpenBitSet(unimpAlign.getSiteCount());
        int swapConflicts=0;
        for (int i = 0; i < unimpAlign.getSiteCount(); i++) {
            if((donorAlign.getMajorAllele(i)!=unimpAlign.getMajorAllele(i))||(donorAlign.getMinorAllele(i)!=unimpAlign.getMinorAllele(i))) {
                swapConflicts++;
                swapMjMnMask.set(i);
            }  
        }
        swapMjMnMask.not();
        System.out.println("swapConflicts"+swapConflicts+" same:"+swapMjMnMask.cardinality());
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(this.unimpAlign);
        Random r=new Random(0);
        int total=0, hybrid=0;
        for (int bt = 0; bt < unimpAlign.getSequenceCount(); bt+=1) {
            int taxaImpCnt=0;
            DonorHypoth[][] regionHypth=new DonorHypoth[blocks][10];
            String name=unimpAlign.getIdGroup().getIdentifier(bt).getFullName();
  //          System.out.printf("Imputing %d:%s ... %n", bt,name);
            long time=System.currentTimeMillis();
            long[] mjT=unimpAlign.getAllelePresenceForAllSites(bt, 0).getBits();
            long[] mnT=unimpAlign.getAllelePresenceForAllSites(bt, 1).getBits();
            for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                int[] resultRange=getBlockWithMinMinorCount(mnT, focusBlock, minMinorCnt);
                total++;
                regionHypth[focusBlock]=getBestInbredDonors(bt, resultRange[0],resultRange[2], focusBlock);
                if(regionHypth[focusBlock][0].getErrorRate()>maximumInbredError) {
                    hybrid++;
                    regionHypth[focusBlock]=getBestDonors(bt, resultRange[0],resultRange[2], focusBlock);
                }
//                regionHypth[focusBlock][0].donor1Taxon=r.nextInt(donorAlign.getSequenceCount());
//                regionHypth[focusBlock][0].donor2Taxon=r.nextInt(donorAlign.getSequenceCount());
                if(resolveMethod==0) setAlignmentWithDonors(regionHypth[focusBlock][0],mna);
            }
            if(resolveMethod==1) solveRegionally(mna, bt, regionHypth);
        }
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile, unImpTargetFile, minMinorCnt));
        if(maskAndTest) s.append(compareSites(mna));
        if(testing>0) s.append(String.format("ParentsRight:%d Wrong %d", parentsRight, parentsWrong));
        System.out.println(s.toString());
        
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
        System.out.println(total+"total : hybrid"+hybrid);
      
       //if we put the share size in the tree map, only remove those at a transiti0n boundary 
       // System.out.printf("L%d R%d p:%g %n",ss.startBlock, ss.endBlock, ss.p);
    }
    
    /**
     * If the target regions has Mendelian errors that it looks for overlapping regional
     * solutions that are better.
     * @param mna
     * @param targetTaxon
     * @param regionHypth 
     */
    private void solveRegionally(MutableNucleotideAlignment mna, int targetTaxon, DonorHypoth[][] regionHypth) {
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            DonorHypoth cbh=regionHypth[focusBlock][0];
//            System.out.printf("%d %d %d %n", regionHypth[focusBlock][0].targetTaxon, 
//                    regionHypth[focusBlock][0].focusBlock, regionHypth[focusBlock][0].mendelianErrors);
            if(regionHypth[focusBlock][0].mendelianErrors==0) {setAlignmentWithDonors(cbh,mna);}
            else {
                int minMendelErrors=cbh.mendelianErrors;
                int currBestBlock=cbh.focusBlock;
                int bestBlockDistance=Integer.MAX_VALUE;
                for (int i = cbh.startBlock; i <= cbh.endBlock; i++) {
                    if((regionHypth[i][0].startBlock<=focusBlock)&&(regionHypth[i][0].endBlock>=focusBlock)) {
                        if(regionHypth[i][0].mendelianErrors<minMendelErrors) {
                            currBestBlock=i;
                            minMendelErrors=regionHypth[i][0].mendelianErrors;
                            bestBlockDistance=Math.abs(i-focusBlock);
                        } 
                        else if((regionHypth[i][0].mendelianErrors==minMendelErrors)&&
                                (Math.abs(i-focusBlock)<bestBlockDistance)) {
                            currBestBlock=i;
                            minMendelErrors=regionHypth[i][0].mendelianErrors;
                            bestBlockDistance=Math.abs(i-focusBlock);
                        }
                    }
                    
                }
//                System.out.printf("c: %d %d %d %n", regionHypth[currBestBlock][0].targetTaxon, 
//                    regionHypth[currBestBlock][0].focusBlock, regionHypth[currBestBlock][0].mendelianErrors);
                setAlignmentWithDonors(regionHypth[currBestBlock][0],mna);
            }
        }
    }
    
    /**
     * Given a start 64 site block, it expands to the left and right until it hits
     * the minimum Minor Site count in the target taxon
     * @param mnT - minor allele bit presence in a series of longs
     * @param focusBlock
     * @param minMinorCnt
     * @return arrays of blocks {startBlock, focusBlock, endBlock}
     */
    private int[] getBlockWithMinMinorCount(long[] mnT, int focusBlock, int minMinorCnt) {
        int minorCnt=Long.bitCount(mnT[focusBlock]);
        int endBlock=focusBlock, startBlock=focusBlock;
        while(minorCnt<minMinorCnt) {
            boolean preferMoveStart=(focusBlock-startBlock<endBlock-focusBlock)?true:false;
            if(startBlock==0) preferMoveStart=false;
            if(endBlock==blocks-1) preferMoveStart=true;
            if(preferMoveStart) {//expand start
                startBlock--;
                minorCnt+=Long.bitCount(mnT[startBlock]);
            } else { //expand end
                endBlock++;
                minorCnt+=Long.bitCount(mnT[endBlock]);  
            } 
        }
        int[] result={startBlock, focusBlock, endBlock};
        return result;
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
    private DonorHypoth[] getBestInbredDonors(int targetTaxon, int startBlock, int endBlock, int focusBlock) {
        long[] mjT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 0).getBits();
        long[] mnT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 1).getBits();
        int[] donors={-1,-1,0};
        double minPropUnmatched=1.0;
        int maxTestSites=0;
        int donorTieCnt=0;
        int[] rDonors=parseDonorsNamesForRegion(unimpAlign.getTaxaName(targetTaxon),startBlock*64);
        if(testing>3) System.out.printf("StartSite %d EndSite %d RealD1 %d RealD2 %d %n",startBlock*64, 
                (endBlock*64+63),rDonors[0],rDonors[1]);
        long[] swapMask=this.swapMjMnMask.getBits();
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        bestDonors.put(1.0, new DonorHypoth());
        
        for (int d1 = 0; d1 < donorAlign.getSequenceCount(); d1++) {
            long[] mj1=donorAlign.getAllelePresenceForAllSites(d1, 0).getBits();
            long[] mn1=donorAlign.getAllelePresenceForAllSites(d1, 1).getBits();
            int mjUnmatched=0;
            int mnUnmatched=0;
            int testSites=0;
            int testTargetMajor=0;
            int testTargetMinor=0;
            for (int i = startBlock; i <= endBlock; i++) {
                long siteMask=swapMask[i]&(mjT[i]|mnT[i])&(mj1[i]|mn1[i]);
                mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i]));
                mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i]));
                testSites+=Long.bitCount(siteMask);
                testTargetMajor+=Long.bitCount(siteMask&mjT[i]);
                testTargetMinor+=Long.bitCount(siteMask&mnT[i]);
            }
            int totalMendelianErrors=mjUnmatched+mnUnmatched;
            double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
            if(testPropUnmatched<bestDonors.lastKey()) {
                DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, startBlock,
                    focusBlock, endBlock, testSites, totalMendelianErrors);
                bestDonors.put(new Double(testPropUnmatched), theDH);
                if(bestDonors.size()>10) bestDonors.remove(bestDonors.lastKey());
            }
            if((testing>1)&&(rDonors[0]==d1)&&(rDonors[1]==d1))
                System.out.printf("Donor %d %d %d %d %d %d %d block %d-%d-%d %n", d1, d1, mjUnmatched, mnUnmatched,
                        testSites, testTargetMajor, testTargetMinor, startBlock, focusBlock, endBlock);

            if((testPropUnmatched<minPropUnmatched)||
                    ((testPropUnmatched==minPropUnmatched)&&(testSites>maxTestSites))) {
                donors[0]=d1;
                donors[1]=d1;
                minPropUnmatched=testPropUnmatched;
                donors[2]=maxTestSites=testSites;
                donorTieCnt=0;
            } else if(testPropUnmatched==minPropUnmatched) donorTieCnt++;
            
        }
        if(testing>1){
            if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {
                System.out.printf("Correct ties:%d bestMatch:%g %n",donorTieCnt,minPropUnmatched);}
            else {System.out.printf("WRONG D1:%d D2:%d ties:%d bestMatch:%g %n",donors[0], donors[1],donorTieCnt,minPropUnmatched);}
        }
        if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {parentsRight++;}
            else {parentsWrong++;}
        DonorHypoth[] result=new DonorHypoth[10];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
        return result;
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
    private DonorHypoth[] getBestDonors(int targetTaxon, int startBlock, int endBlock, int focusBlock) {
        long[] mjT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 0).getBits();
        long[] mnT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 1).getBits();
        int[] donors={-1,-1,0};
        double minPropUnmatched=1.0;
        int maxTestSites=0;
        int donorTieCnt=0;
        int[] rDonors=parseDonorsNamesForRegion(unimpAlign.getTaxaName(targetTaxon),startBlock*64);
        if(testing>3) System.out.printf("StartSite %d EndSite %d RealD1 %d RealD2 %d %n",startBlock*64, 
                (endBlock*64+63),rDonors[0],rDonors[1]);
        long[] swapMask=this.swapMjMnMask.getBits();
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        bestDonors.put(1.0, new DonorHypoth());
        
        for (int d1 = 0; d1 < donorAlign.getSequenceCount(); d1++) {
            long[] mj1=donorAlign.getAllelePresenceForAllSites(d1, 0).getBits();
            long[] mn1=donorAlign.getAllelePresenceForAllSites(d1, 1).getBits();
            for (int d2 = d1; d2 < donorAlign.getSequenceCount(); d2++) {
                long[] mj2=donorAlign.getAllelePresenceForAllSites(d2, 0).getBits();
                long[] mn2=donorAlign.getAllelePresenceForAllSites(d2, 1).getBits();
                int mjUnmatched=0;
                int mnUnmatched=0;
                int testSites=0;
                int testTargetMajor=0;
                int testTargetMinor=0;
                for (int i = startBlock; i <= endBlock; i++) {
                    long siteMask=swapMask[i]&(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
                    mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
                    mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
                    testSites+=Long.bitCount(siteMask);
                    testTargetMajor+=Long.bitCount(siteMask&mjT[i]);
                    testTargetMinor+=Long.bitCount(siteMask&mnT[i]);
                }
                int totalMendelianErrors=mjUnmatched+mnUnmatched;
                double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
                if(testPropUnmatched<bestDonors.lastKey()) {
                    DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, startBlock,
                        focusBlock, endBlock, testSites, totalMendelianErrors);
                    bestDonors.put(new Double(testPropUnmatched), theDH);
                    if(bestDonors.size()>10) bestDonors.remove(bestDonors.lastKey());
                }
                if((testing>1)&&(rDonors[0]==d1)&&(rDonors[1]==d2))
                    System.out.printf("Donor %d %d %d %d %d %d %d block %d-%d-%d %n", d1, d2, mjUnmatched, mnUnmatched,
                            testSites, testTargetMajor, testTargetMinor, startBlock, focusBlock, endBlock);
                
                if((testPropUnmatched<minPropUnmatched)||
                        ((testPropUnmatched==minPropUnmatched)&&(testSites>maxTestSites))) {
                    donors[0]=d1;
                    donors[1]=d2;
                    minPropUnmatched=testPropUnmatched;
                    donors[2]=maxTestSites=testSites;
                    donorTieCnt=0;
                } else if(testPropUnmatched==minPropUnmatched) donorTieCnt++;
                
            }
            
        }
        if(testing>1){
            if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {
                System.out.printf("Correct ties:%d bestMatch:%g %n",donorTieCnt,minPropUnmatched);}
            else {System.out.printf("WRONG D1:%d D2:%d ties:%d bestMatch:%g %n",donors[0], donors[1],donorTieCnt,minPropUnmatched);}
        }
        if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {parentsRight++;}
            else {parentsWrong++;}
        DonorHypoth[] result=new DonorHypoth[10];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
        return result;
    }
  
    
    /**
     * Takes a donor hypothesis and applies it to the output alignment 
     * @param theDH
     * @param mna 
     */
    private void setAlignmentWithDonors(DonorHypoth theDH, MutableNucleotideAlignment mna) {
        int startSite=theDH.focusBlock*64;
        int endSite=startSite+63;
        if(endSite>=unimpAlign.getSiteCount()) endSite=unimpAlign.getSiteCount()-1;
        for(int cs=startSite; cs<=endSite; cs++) {
            byte bD1=donorAlign.getBase(theDH.donor1Taxon, cs);
            byte bD2=donorAlign.getBase(theDH.donor2Taxon, cs);
            byte donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE;
            if(bD1==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD2;}
            else if(bD2==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD1;}
            else {donorEst=(byte)((bD1&highMask)|(bD2&lowMask));
            }
            //need to check whether the heterozygote is put together properly
            //need to change to 
            mna.setBase(theDH.targetTaxon, cs, donorEst);
//            if(mna.getBase(theDH.targetTaxon, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                    mna.setBase(theDH.targetTaxon, cs, donorEst);
////                    impSiteCnt++;
////                    taxaImpCnt++;
//                }
        } //end of cs loop

    }
    
    private int[] parseDonorsNamesForRegion(String taxonName, int site) {
        String[] sb=taxonName.split("\\|");
        int[] donors={-1,-1};
        for (int i = 1; i < sb.length; i++) {
            String[] ss=sb[i].split("s");
            int leftSite=Integer.parseInt(ss[1]);
            if(leftSite<=site) {
                String[] sd=ss[0].split("_");
                donors[0]=Integer.parseInt(sd[0]);
                donors[1]=Integer.parseInt(sd[1]);
              //  break;
            }
        }
        return donors;
    }
    
    /**
     * Used for testing imputation accuracy.  Sets every sampIntensity sites to missing
     * so they can be compared at the end of imputation.
     * @param sampIntensity 
     */
    private void maskSites(int sampIntensity) {
        System.out.println("Beginning to mask sites");
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(unimpAlign);
        int maxSites=mna.getSequenceCount()*((mna.getSiteCount()/sampIntensity)+1);
        hSite=new int[maxSites];
        hTaxon=new int[maxSites];
        hState=new byte[maxSites];
        int cnt=0;
        for (int t = 0; t < mna.getSequenceCount(); t++) {
            for (int s = t%sampIntensity; s < mna.getSiteCount(); s+=sampIntensity) {
                hSite[cnt]=s;
                hTaxon[cnt]=t;
                hState[cnt]=mna.getBase(t, s);
                mna.setBase(t, s, Alignment.UNKNOWN_DIPLOID_ALLELE);
                cnt++;
            }
 //           System.out.println(t+":"+cnt);
        }
        maskSitCnt=cnt;
        mna.clean();
        compareSites(unimpAlign);
        unimpAlign=BitAlignment.getInstance(mna, false);
        compareSites(unimpAlign);
        System.out.println("Sites masked");
    }
    
    private String compareSites(Alignment a) {
        int missingCnt=0, correctCnt=0, errorCnt=0, notImp=0, hetCnt=0;
        for (int i = 0; i < maskSitCnt; i++) {
            if(hState[i]==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                missingCnt++;
                continue;
            }
            byte impb=a.getBase(hTaxon[i], hSite[i]);
            if(AlignmentUtils.isHeterozygous(impb)) {
                hetCnt++;
                byte[] orig=AlignmentUtils.getDiploidValues(hState[i]);
                byte[] imphet=AlignmentUtils.getDiploidValues(impb);
                if((orig[0]==imphet[0])||(orig[0]==imphet[1])||(impb==hState[i])) {correctCnt++;}
                else {errorCnt++;}
            } else if(impb==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                notImp++;
            } else if(impb==hState[i]) {
                correctCnt++;
            } else {errorCnt++;}
        }
        double errRate=(double)errorCnt/(double)(errorCnt+correctCnt);
        return String.format("Missing: %d Het: %d NotImp: %d Error: %d Correct: %d ErrorRate: %g ", 
                missingCnt, hetCnt, notImp, errorCnt, correctCnt, errRate);
        
    }
    
    private String reportTaxaMakeUp(ArrayList<Integer>[] data) {
        StringBuilder s=new StringBuilder();
        for (int i = 0; i < data[0].size(); i++) {
            s.append(data[0].get(i));
            s.append(":");
            s.append(data[1].get(i));
            s.append("\t");
        }
        return s.toString();
    }
   
    
    /**
     * 
     * @param iMj
     * @param iMn
     * @param jMj
     * @param jMn
     * @return [0]= number of sites in comparison,
     * [1]=number of minor alleles in comparison from taxon i  
     * [2]=sites that agree
     */
    private int[] getIdentity(long iMj, long iMn, long jMj, long jMn) {
        int[] results=new int[3];
        long iMnjMn=iMn&jMn;
        long iMnjMj=iMn&jMj;
        long iMnComps=iMnjMn|iMnjMj;
        long sameL=(iMj&jMj)|(iMnjMn);
        long diffL=(iMj&jMn)|(iMnjMj);
        long hetsL=sameL&diffL;
        int same=(int)BitUtil.pop(sameL);
        int diff=(int)BitUtil.pop(diffL);
        int hets=(int)BitUtil.pop(hetsL);
        results[1]=(int)BitUtil.pop(iMnComps);  
        int sum=same+diff+hets;
        results[0]=sum-(2*hets);
        results[2]=same-(hets/2); //check the minus sign
        return results;
    }
    
    public static void createSynthetic(String donorFile, String unImpTargetFile, int blockSize,
            double propPresent, double inbreedingF, int taxaNumber) {
        Alignment a=ImportUtils.readFromHapmap(donorFile, (ProgressListener)null);
        System.out.printf("Read %s Sites %d Taxa %d %n", donorFile, a.getSiteCount(), a.getSequenceCount());
        MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a, taxaNumber, a.getSiteCount());
        Random r=new Random();
        for (int t = 0; t < taxaNumber; t++) {
            StringBuilder tName=new StringBuilder("ZM"+t);
            int p1=r.nextInt(a.getSequenceCount());
            int p2=r.nextInt(a.getSequenceCount());
            for (int b = 0; b < a.getSiteCount(); b+=blockSize) {  //change to bp?
                if(r.nextDouble()<0.5) {p1=r.nextInt(a.getSequenceCount());}
                else {p2=r.nextInt(a.getSequenceCount());}
                if(r.nextDouble()<inbreedingF) {p2=p1;}  //inbreeding
                if(p2<p1) {int temp=p1; p1=p2; p2=temp;}
                tName.append("|"+p1+"_"+p2+"s"+b);
                for (int s = b; (s < b+blockSize) && (s<a.getSiteCount()); s++) {
                    double presentRoll=r.nextDouble();
                    byte p1b=a.getBase(p1, s);
                    if(AlignmentUtils.isHeterozygous(p1b)==true) {p1b=Alignment.UNKNOWN_DIPLOID_ALLELE;}
                    byte p2b=a.getBase(p2, s);
                    if(AlignmentUtils.isHeterozygous(p2b)==true) {p2b=Alignment.UNKNOWN_DIPLOID_ALLELE;}
                    if(presentRoll<(propPresent*propPresent)) {
                        if(p1b==Alignment.UNKNOWN_DIPLOID_ALLELE) {mna.setBase(t, s, p2b);}
                        else if(p2b==Alignment.UNKNOWN_DIPLOID_ALLELE) {mna.setBase(t, s, p1b);}
                        else {
                            byte het=AlignmentUtils.getUnphasedDiploidValueNoHets(p1b, p2b);
                            mna.setBase(t, s, het);
                        }
                    } else if(presentRoll<propPresent) {
                        if(r.nextDouble()<0.5) {
                            mna.setBase(t, s, p1b);
                        } else {
                            mna.setBase(t, s, p2b);
                        }
                        //if(AlignmentUtils.isHeterozygous(mna.getBase(t, s))==true) {mna.setBase(t, s, Alignment.UNKNOWN_DIPLOID_ALLELE);}
                    } else {
                        mna.setBase(t, s, Alignment.UNKNOWN_DIPLOID_ALLELE);
                    }
                    if(r.nextDouble()<0.000) {  //method to add error to simulation 0.002 is reasonable
                        mna.setBase(t, s, AlignmentUtils.getUnphasedDiploidValueNoHets(a.getMinorAllele(s),a.getMinorAllele(s)));
                    }
                    if((p1==p2)&&(AlignmentUtils.isHeterozygous(mna.getBase(t, s))==true)) {
                        System.out.printf("%d %d %d %d %d %d %s %s %s %n", p1, p2, p1b, p2b, a.getBase(p1, s), a.getBase(p2, s), a.getBaseAsString(p1, s), a.getBaseAsString(p2, s), mna.getBaseAsString(t, s));
                    }
                }//end of site        
            } //end of blocks
            System.out.println(tName.toString());
            if(mna.getSequenceCount()==t) {mna.addTaxon(new Identifier(tName.toString()));}
            else {mna.setTaxonName(t, new Identifier(tName.toString()));}
        }
        mna.clean();
        ExportUtils.writeToHapmap(mna, false, unImpTargetFile, '\t', null);
    }
    
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
      String root="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/";
//        String root="/Volumes/LaCie/build20120110/imp/";

        String donorFile=root+"NAMfounder20120110.imp.hmp.txt";
 //       String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
        String unImpTargetFile=root+"ZeaSyn20120110.hmp.txt";
        String impTargetFile=root+"ZeaSyn20120110.imp.hmp.txt";

        boolean buildInput=true;
        boolean filterTrue=true;
        if(buildInput) {createSynthetic(donorFile, unImpTargetFile, 2000, 0.4, 0.5, 1000);}

        KnownParentMinorWindowImputation e64NNI=new KnownParentMinorWindowImputation(donorFile,
                unImpTargetFile, impTargetFile,15,1);


        
//        for (int recSize = 512; recSize < 4000; recSize+=(recSize/2)) {
//            for (int mm = 5; mm < 40; mm+=5) {
//                System.out.println("Rec size"+recSize);
//                unImpTargetFile=root+recSize+"ZeaSyn20120110.hmp.txt";
//                if(buildInput) {createSynthetic(donorFile, unImpTargetFile, recSize, 0.4, -1, 1000);}
//                e64NNI=new KnownParentMinorWindowImputation(donorFile, unImpTargetFile, impTargetFile,mm,1);
//                }
//        }
    }
    
}
