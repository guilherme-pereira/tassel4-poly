/*
 * Expanding64NNImputation
 */
package net.maizegenetics.pal.popgen;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Random;
import java.util.TreeMap;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.statistics.ApproxFastChiSquareDistribution;
import net.maizegenetics.pipeline.EdTests;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ProgressListener;

/**
 * Finds the nearest neighbor for every 64 site window.  In case of ties, it 
 * extends to neighboring 64bp windows.  
 * 
 * @author edbuckler
 */
public class Expanding64NNImputation {
    private Alignment ldAlign;
    int minSites=256;
    int maxWindow=2048/64;
    double minIdentityDiff=0.01;
    int[][][] null64Share; //region, site cnt, siteIdentity
    float[][][] null64ShareProb; //region, site cnt, siteIdentity
    int blocks=-1;
    int[] hSite, hTaxon;
    byte[] hState;
    int maskSitCnt=0;
    int maxNN=50;
    double minProb=0.0001;
    boolean maskAndTest=true;
    ApproxFastChiSquareDistribution fcs=new ApproxFastChiSquareDistribution(1000,200);

    public Expanding64NNImputation(Alignment inldAlign, String exportFile) {
        this.ldAlign = AlignmentUtils.optimizeForTaxa(inldAlign);
        if(maskAndTest) maskSites(300);
        blocks=ldAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
//        this.createNull64Share(ldAlign, 500);
        this.createNull64Share(ldAlign);
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(this.ldAlign);
        int impSiteCnt=0;
        for (int bt = 0; bt < ldAlign.getSequenceCount(); bt++) {
            int taxaImpCnt=0;
            String name=ldAlign.getIdGroup().getIdentifier(bt).getFullName();
            System.out.printf("Imputing %d:%s ...", bt,name);
            long time=System.currentTimeMillis();
            float[][] idp=getTaxaIdentityProbMatrix(bt);
            System.out.printf("IDmatrixTime %d ", System.currentTimeMillis()-time);
            for (int x = 0; x < idp[0].length; x++) {
                int startSite=x*64;
                int endSite=startSite+63;
                if(endSite>=ldAlign.getSiteCount()) endSite=ldAlign.getSiteCount()-1;
                TreeMap<Double, ShareSize> bestTaxa=new TreeMap<Double, ShareSize>();
                for (int t = 0; t < idp.length; t++) {
                    ShareSize xss=new ShareSize(bt,t, x, x);
                    ShareSize fss=getMaxShare(idp,xss);
                    if(fss.p>minProb) continue;
                    if((bestTaxa.size()<maxNN)||(fss.p<bestTaxa.lastEntry().getKey())) {
                        bestTaxa.put(fss.p, fss);
                        if(bestTaxa.size()>maxNN) bestTaxa.remove(bestTaxa.lastEntry().getKey());
                    }
                    //System.out.printf("%g\t",idp[t][x]);   
                }
//                byte[] calls=this.consensusCallBit(bt, x, bestTaxa, false, 0.66, 4, true, false);
//                for(int cs=startSite; cs<=endSite; cs++) {
//                        if(mna.getBase(bt, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                            byte nb=calls[cs-startSite];
//                            if(nb!=Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                                mna.setBase(bt, cs, nb);
//                                impSiteCnt++;
//                                taxaImpCnt++;
//                            }
//                        }
//                    }
//                for(ShareSize c: bestTaxa.values()) {
//                    int ct=c.compTaxon;
//                    for(int cs=startSite; cs<=endSite; cs++) {
//                        if(mna.getBase(bt, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                            byte nb=ldAlign.getBase(ct, cs);
//                            if(nb!=Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                                mna.setBase(bt, cs, nb);
//                                impSiteCnt++;
//                                taxaImpCnt++;
//                            }
//                        }
//                    }
//                }
                
                for(int cs=startSite; cs<=endSite; cs++) {
                    if(mna.getBase(bt, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                        mna.setBase(bt, cs, getBestBase(mna,bestTaxa.values(),cs));
                        impSiteCnt++;
                        taxaImpCnt++;
                    }
                }
            }
            System.out.printf("Finished %d Imp %d %d %n", System.currentTimeMillis()-time, impSiteCnt, taxaImpCnt);
            if(bt%10==0) compareSites(mna);
        }
        if(maskAndTest) compareSites(mna);
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
      
       //if we put the share size in the tree map, only remove those at a transiti0n boundary 
       // System.out.printf("L%d R%d p:%g %n",ss.left, ss.right, ss.p);
    }
    
    private byte getBestBase(Alignment mna, Collection<ShareSize> bestTaxa, int cs) {
        byte mjA=mna.getMajorAllele(cs);
        mjA=(byte)((mjA<<4)|mjA);
        byte mnA=mna.getMinorAllele(cs);
        mnA=(byte)((mnA<<4)|mnA);
        int mjCnt=0, mnCnt=0, unkCnt=0;
        double mjLnSum=0, mnLnSum=0, unkLnSum=0;        
        for(ShareSize c: bestTaxa) {
            int ct=c.compTaxon;
            byte nb=ldAlign.getBase(ct, cs);
            if(nb==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                unkCnt++;
          //      unkLnSum+=Math.log(c.p);
            } else if(nb==mjA) {
                mjCnt++;
          //      mjLnSum+=Math.log(c.p);
            } else if(nb==mnA) {
                mnCnt++;
           //     mnLnSum+=Math.log(c.p);
            }
        }
        if((mnCnt>(mjCnt+1))&&(mnCnt>0)) return mnA;
        if(((1+mnCnt)<mjCnt)&&(mjCnt>0)) return mjA;
        return Alignment.UNKNOWN_DIPLOID_ALLELE;
    }
    
    private byte[] consensusCallBit(int taxon, int block, TreeMap<Double,ShareSize> taxa,
            boolean callhets, double majority, int minCount, boolean ignoreKnownBases, boolean imputeGaps) {
        int[] taxaIndex=new int[taxa.size()];
        ArrayList<ShareSize> taxaList=new ArrayList(taxa.values());
        for (int t = 0; t < taxaIndex.length; t++) {
            taxaIndex[t]=taxaList.get(t).compTaxon;
        }
        short[][] siteCnt=new short[2][64];
//        double[] sumExpPresent=new double[endBase-startBase];
//        int[] sumNNxSitePresent=new int[endBase-startBase];
//        int sumNNPresent=0;
     //   for (int t = 0; t < taxaIndex.length; t++) sumNNPresent+=presentCntForTaxa[taxaIndex[t]];

            int currWord=block;
            int callSite=0;
            for (int t = 0; t < taxaIndex.length; t++) {
//                long bmj=ldAlign.getAllelePresenceForSitesBlock(taxaIndex[t], 0,currWord, currWord+1)[0];
//                long bmn=ldAlign.getAllelePresenceForSitesBlock(taxaIndex[t], 1,currWord, currWord+1)[0];
                long bmj=ldAlign.getAllelePresenceForAllSites(taxaIndex[t], 0).getBits()[block];
                long bmn=ldAlign.getAllelePresenceForAllSites(taxaIndex[t], 1).getBits()[block];
                int cs=callSite;
                for (int j = 0; j < 64; j++) {
                    boolean presentFlag=false;
                    if((bmj & 0x01)!=0) {siteCnt[0][cs]++; presentFlag=true;}
                    bmj=bmj>>1;
                    if((bmn & 0x01)!=0) {siteCnt[1][cs]++; presentFlag=true;}
                    bmn=bmn>>1;
       //             sumExpPresent[cs]+=presentProp[taxaIndex[t]];
//                    if(presentFlag) sumNNxSitePresent[cs]++;
                    cs++;          
                }
            }
 //       System.out.println("Bit:"+Arrays.toString(siteCnt[0]));
        byte[] calls=new byte[64];
        Arrays.fill(calls, Alignment.UNKNOWN_DIPLOID_ALLELE);
        int startSite=block*64;
        int endSite=startSite+63;
        for (int alignS = startSite; alignS <= endSite; alignS++) {
            int callS=alignS-startSite;
            byte ob=ldAlign.getBase(taxon,alignS);
            
            if(ignoreKnownBases) ob=Alignment.UNKNOWN_DIPLOID_ALLELE;
            calls[callS]=ob;
            byte mj=ldAlign.getMajorAllele(alignS);
            byte mn=ldAlign.getMinorAllele(alignS);
            mj=(byte)((mj<<4)+mj);
            mn=(byte)((mn<<4)+mn);
            
            int totalCnt=siteCnt[0][callS]+siteCnt[1][callS];
 //           double expPres=sumExpPresent[callS]/(double)taxaIndex.length;
            
            if(totalCnt<minCount) continue;  //no data leave missing
            if((double)siteCnt[0][callS]/(double)totalCnt>majority) {
                if((ob!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(ob!=mj)) {calls[callS]=Alignment.UNKNOWN_DIPLOID_ALLELE;}
                else {calls[callS] = mj;}
            }
            else if((double)siteCnt[1][callS]/(double)totalCnt>majority) {
                if((ob!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(ob!=mn)) {calls[callS]=Alignment.UNKNOWN_DIPLOID_ALLELE;}
                else {calls[callS] = mn;}
            }
            else if(callhets) {
//                byte[] snpValue={mj,mn};
//                byte het=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(snpValue);
//                calls[callS]=het;
            }
 //           System.out.printf("Taxon:%d orig:%d mj:%d mn:%d call:%d %s %n", taxon, ldAlign.getBase(taxon,alignS),  mj, mn, calls[callS], AlignmentUtils.isHeterozygous(calls[callS]));
        }
        return calls;
    }
    
    private void maskSites(int sampIntensity) {
        System.out.println("Beginning to mask sites");
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(ldAlign);
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
        compareSites(ldAlign);
        ldAlign=BitAlignment.getInstance(mna, false);
        compareSites(ldAlign);
        System.out.println("Sites masked");
    }
    
    private void compareSites(Alignment a) {
        int missingCnt=0, correctCnt=0, errorCnt=0, notImp=0, hetCnt=0;
        for (int i = 0; i < maskSitCnt; i++) {
            if(hState[i]==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                missingCnt++;
                continue;
            }
            byte impb=a.getBase(hTaxon[i], hSite[i]);
            if(AlignmentUtils.isHeterozygous(impb)) {
                hetCnt++;
            } else if(impb==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                notImp++;
            } else if(impb==hState[i]) {
                correctCnt++;
            } else {errorCnt++;}
        }
        double errRate=(double)errorCnt/(double)(errorCnt+correctCnt);
        System.out.printf("Missing: %d Het: %d NotImp: %d Error: %d Correct: %d ErrorRate: %g %n", 
                missingCnt, hetCnt, notImp, errorCnt, correctCnt, errRate);
        
    }
    
    private ShareSize getMaxShare(float[][] idp, ShareSize currShare) {
        if(currShare.left==currShare.right) {
            currShare.fsum=idp[currShare.compTaxon][currShare.left];
     //       currShare.p=1.0-ChiSquareDistribution.cdf(currShare.fsum, 2);
            currShare.p=1.0-fcs.cdfFastApprox(currShare.fsum, 2);
        }
        double tL=-1, tR=-1;
        if(currShare.left>0) {
            tL=idp[currShare.compTaxon][currShare.left-1];
        }
        if(currShare.right<idp[0].length-1) {
            tR=idp[currShare.compTaxon][currShare.right+1];
        }
        if(tL>tR) {
            double testFsum=currShare.fsum+tL;
         //   double testp=1.0-ChiSquareDistribution.cdf(testFsum, currShare.df+2);
            double testp=1.0-fcs.cdfFastApprox(testFsum, currShare.df+2);
            if(testp<currShare.p) {
                currShare.moveLeft(testp, testFsum);
                return getMaxShare(idp, currShare);
            }
        } else {
            double testFsum=currShare.fsum+tR;
//            double testp=1.0-ChiSquareDistribution.cdf(testFsum, currShare.df+2);
            double testp=1.0-fcs.cdfFastApprox(testFsum, currShare.df+2);
            if(testp<currShare.p) {
                currShare.moveRight(testp, testFsum); 
                return getMaxShare(idp, currShare);
            }
        }
        return currShare;
    }

    private void createNull64Share(Alignment a, int maxSampling) {
        a = AlignmentUtils.optimizeForTaxa(a);
        System.out.println("Creating the null distribution");
        IBSDistanceMatrix dm=new IBSDistanceMatrix(a,100,null);
        System.out.printf("Distances estimated. Mean:%g %n", dm.meanDistance());
        double meanDist=dm.meanDistance();
        null64Share=new int[blocks][65][65];
        for (int i = 0; i < blocks; i++) {
            for (int j = 0; j < null64Share[0].length; j++) {
                for (int k = 0; k <=j; k++) {
                    null64Share[i][j][k]=1;
                }                
            } 
        }
        Random r=new Random(0);
      //  int samplingPerTaxon=maxSampling/(a.getSequenceCount()*a.getSequenceCount()/2);
        System.out.println("samplingPerTaxonContrast: "+maxSampling);
        for (int t1 = 0; t1 < a.getSequenceCount(); t1++) {
            long[] iMj=a.getAllelePresenceForAllSites(t1, 0).getBits();
            long[] iMn=a.getAllelePresenceForAllSites(t1, 1).getBits();
            for (int samp = 0; samp < maxSampling; samp++) {
                int d=0;
                int t2=r.nextInt(a.getSequenceCount());
                while(dm.getDistance(t1, t2)<meanDist) {
                    t2=r.nextInt(a.getSequenceCount());
                }
 //               int t2=r.nextInt(a.getSequenceCount());
 //               int t2=getIndexOfMaxDistance(dm, t1);
                if(t1==t2) continue;
                long[] jMj=a.getAllelePresenceForAllSites(t2, 0).getBits();
                long[] jMn=a.getAllelePresenceForAllSites(t2, 1).getBits();
                for (int sN = 0; sN < blocks; sN++) {
               //     int br=r.nextInt(numBins);
                    int b=sN;
                    int[] results=this.getIdentity(iMj[b], iMn[b], jMj[b], jMn[b]);
              //      if((sN==0)&&((t1==39)||(t2==39))) System.out.printf("%d %d %d %s %n",sN, t1, t2, Arrays.toString(results));
                    null64Share[sN][results[0]][results[2]]++;
                }
            }
        }
        null64ShareProb=new float[blocks][65][65];
        for (int i = 0; i < blocks; i++) {
            for (int j = 0; j < null64Share[0].length; j++) {
                int sum=0, bsum=0;
                for (int k = 0; k <=j; k++) {
                    sum+=null64Share[i][j][k];
                }
                for (int k = j; k >=0; k--) {
                    bsum+=null64Share[i][j][k];
                  //  null64ShareProb[i][j][k]=(float)bsum/(float)sum;
                    null64ShareProb[i][j][k]=(float)(-2*Math.log((double)bsum/(double)sum));
                }  
 //               System.out.printf("%d %g %d %s %n", i, 0.5, j, Arrays.toString(null64ShareProb[i][j]));
            } 
        }
    }
    
    private void createNull64Share(Alignment a) {
        a = AlignmentUtils.optimizeForTaxa(a);
        System.out.println("Creating the SBitAlignment distribution");
        Alignment sbit=BitAlignment.getInstance(a, true);
        System.out.printf("SBitAlignment created %n");

        null64ShareProb=new float[blocks][65][66];
        net.maizegenetics.pal.math.Binomial bn=new net.maizegenetics.pal.math.Binomial();
        for (int i = 0; i < blocks; i++) {
            double mafSum=0;
            int cnt=0;
            for (int s = i*64; (s < (i+1)*64)&&(s<sbit.getSiteCount()); s++) { 
                mafSum+=sbit.getMinorAlleleFrequency(s);
                cnt++;
            }
            double avgMAF=mafSum/(double)cnt;
            for (int j = 1; j < null64ShareProb[0].length; j++) {
                Binomial binomFunc=new Binomial(j, 1.0-avgMAF, new RandomJava());
                Arrays.fill(null64ShareProb[i][j], 0);
                for (int k = j; k >=0; k--) {
                    null64ShareProb[i][j][k]=(float)(binomFunc.pdf(k))+null64ShareProb[i][j][k+1];
                }
                for (int k = j; k >=0; k--) {
                    if(null64ShareProb[i][j][k]>1) null64ShareProb[i][j][k]=1;
                    null64ShareProb[i][j][k]=-2*(float)Math.log(null64ShareProb[i][j][k]);
                }
  //              System.out.printf("%d %g %d %s %n", i, avgMAF, j, Arrays.toString(null64ShareProb[i][j]));
            } 
        }
    }
    
    private int getIndexOfMaxDistance(DistanceMatrix dm, int compTaxon) {
        int resultTaxon=compTaxon;
        double maxDist=dm.getDistance(compTaxon, compTaxon);
        for (int i = 0; i < dm.getSize(); i++) {
            if(dm.getDistance(compTaxon, i)>maxDist) {
                maxDist=dm.getDistance(compTaxon, i);
                resultTaxon=i;
            }   
        }
        return resultTaxon;
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
    
    private float[][] getTaxaIdentityProbMatrix(int taxa) {
        long[] iMj=ldAlign.getAllelePresenceForAllSites(taxa, 0).getBits();
        long[] iMn=ldAlign.getAllelePresenceForAllSites(taxa, 1).getBits();
        int sections=iMj.length;
        float[][] result=new float[ldAlign.getSequenceCount()][sections];
        for (int t = 0; t < ldAlign.getSequenceCount(); t++) {
            long[] jMj=ldAlign.getAllelePresenceForAllSites(t, 0).getBits();
            long[] jMn=ldAlign.getAllelePresenceForAllSites(t, 1).getBits();
            for(int x=0; x<sections; x++) {
                int[] results=this.getIdentity(iMj[x], iMn[x], jMj[x], jMn[x]);
                result[t][x]=null64ShareProb[x][results[0]][results[2]];
            }
        }
        return result;
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
    
    
    public static void main(String[] args) {
//        String root="/Users/edbuckler/SolexaAnal/bigprojection/";
        String root="/Volumes/LaCie/bigprojection/";
//        String gFile=root+"282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c10.hmp.txt";
//        String mFile=root+"282merge.chr10.hmp.txt";
//        String exFile=root+"282merge_01.chr10.imp.hmp.txt";
        String gFile=root+"Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c10.hmp.txt";
        String mFile=root+"10pctZeamerge.chr10.hmp.txt";
        String exFile=root+"10pctZeamerge_01.chr10.imp.hmp.txt";
        
        String hFile=root+"maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt";

        boolean buildInput=false;
        boolean filterTrue=true;
        
        if(buildInput) {
            Alignment a=ImportUtils.readFromHapmap(gFile, (ProgressListener)null);
            System.out.println("GBS Map Read");
            if(filterTrue) a=FilterAlignment.getInstance(a, 0, a.getSiteCount()/10);
            Alignment gbsMap=BitAlignment.getInstance(a, false);
            
  //          TBitAlignment gbsMap=TBitAlignment.getInstance(ImportUtils.readFromHapmap(gFile, (ProgressListener)null));
            System.out.println("GBS converted and filtered");
    //        SBitAlignment hapMap=(SBitAlignment)readGZOfSBit(hapFileAGP1, true);
            Alignment hapMap=ImportUtils.readFromHapmap(hFile, true, (ProgressListener)null);
            System.out.println("HapMap Read");
            hapMap=EdTests.fixHapMapNames(hapMap);  //adds tags so that HapMapNames are recognizable
            System.out.println("HapMap Names Fixed");
            MutableNucleotideAlignment mna=EdTests.combineAlignments(hapMap, gbsMap);
            System.out.println("HapMap and GBS combined");
            mna.clean();
            ExportUtils.writeToHapmap(mna, false, mFile, '\t', null);
        }
        Alignment mergeMap=ImportUtils.readFromHapmap(mFile, false, (ProgressListener)null);
        Expanding64NNImputation e64NNI=new Expanding64NNImputation(mergeMap, exFile);
    //    TBitAlignment mergeMap=TBitAlignment.getInstance(mna);
    }
    
}

class ShareSize {
    int baseTaxon=-1;
    int compTaxon=-1;
    int left=-1;
    int right=-1;
    double p=1;
    double fsum=0;
    int df=0;

    public ShareSize(int baseTaxon, int compTaxon, int left, int right, double p, double fsum) {
        this(baseTaxon, compTaxon, left, right);
        this.p=p;
        this.fsum=fsum;
    }
    
    public ShareSize(int baseTaxon, int compTaxon, int left, int right) {
        this.baseTaxon=baseTaxon;
        this.compTaxon=compTaxon;
        this.left=left;
        this.right=right;
        df=(right-left+1)*2;
    }
    
    public void moveLeft(double p, double fsum) {
        left--;
        this.p=p;
        this.fsum=fsum;
        df=df+2;
    }
    
    public void moveRight(double p, double fsum) {
        right++;
        this.p=p;
        this.fsum=fsum;
        df=df+2;
    }
    
    public String toString() {
        return String.format("BTx:%d CTx:%d L:%d R:%d FSum:%g P:%g ", baseTaxon,
                compTaxon, left, right, fsum, p);
    }
}
