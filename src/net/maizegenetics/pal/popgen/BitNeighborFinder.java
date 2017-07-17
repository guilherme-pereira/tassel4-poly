/*
 * BitNeighborFinder
 */
package net.maizegenetics.pal.popgen;

import java.util.ArrayList;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ProjectionAlignment;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitUtil;

/**
 * Finds the nearest neighbor for every 64 site window.  In case of ties, it 
 * extends to neighboring 64bp windows.  It can be restricted to only find 
 * neighbors for a taxa within a high density genotype set.  This is useful for 
 * projection.  This class is expected to be used with Projection alignment.
 * 
 * @author edbuckler
 */
public class BitNeighborFinder {
    private byte[][] same, diff, hets;
    IdGroup hdTargetID;
    int[] hdTarget;
    private Alignment ldAlign;
    int minSites=256;
    int maxWindow=2048/64;
    double minIdentityDiff=0.01;
    ProjectionAlignment pa;
    

    public BitNeighborFinder(IdGroup hdTargetID, Alignment ldAlign, Alignment hdAlign) {
        this.hdTargetID = hdTargetID;
        this.ldAlign = AlignmentUtils.optimizeForTaxa(ldAlign);
        initHDTargets();
        pa=new ProjectionAlignment(hdAlign,ldAlign.getIdGroup());
        long prevTime = System.currentTimeMillis();
        for (int i = 0; i < ldAlign.getSequenceCount(); i++) {
            ArrayList<Integer>[] result=findHDTaxa(i);
            int[] posBreaks=new int[result[0].size()];
            int[] hdTaxa=new int[result[0].size()];
            for (int b = 0; b < hdTaxa.length; b++) {
                posBreaks[b]=result[0].get(b);
                hdTaxa[b]=result[1].get(b);
            }
           // pa.setCompositionOfTaxon(i, posBreaks, hdTaxa);
         //   pa.setCompositionOfTaxon(ldAlign.getIdGroup().getIdentifier(i).getFullName(), posBreaks, hdTaxa);
            System.out.printf("%s\t%d\t%s %n",ldAlign.getTaxaName(i),result[0].size(),reportTaxaMakeUp(result));
        }
        long currentTime = System.currentTimeMillis();
        System.out.println("Time to load alleles: " + ((currentTime - prevTime) / 1000));
        
    }

    public ProjectionAlignment getPa() {
        return pa;
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
    
    private ArrayList<Integer>[] findHDTaxa(int taxa) {
 //       System.out.print(ldAlign.getTaxaName(taxa)+":");
        ArrayList<Integer>  posBreaks=new ArrayList();
        ArrayList<Integer>  hdTaxa=new ArrayList();
        long[] iMj=ldAlign.getAllelePresenceForAllSites(taxa, 0).getBits();
        long[] iMn=ldAlign.getAllelePresenceForAllSites(taxa, 1).getBits();
        int sections=iMj.length;
        same=new byte[hdTarget.length][sections];
        diff=new byte[hdTarget.length][sections]; 
        hets=new byte[hdTarget.length][sections];
        for (int t = 0; t < hdTarget.length; t++) {
            if(hdTarget[t]<0) continue;
            long[] jMj=ldAlign.getAllelePresenceForAllSites(hdTarget[t], 0).getBits();
            long[] jMn=ldAlign.getAllelePresenceForAllSites(hdTarget[t], 1).getBits();
            int sameCnt=0, diffCnt=0, hetCnt=0;
            for(int x=0; x<sections; x++) {
                long sameL=(iMj[x]&jMj[x])|(iMn[x]&jMn[x]);
                long diffL=(iMj[x]&jMn[x])|(iMn[x]&jMj[x]);
                long hetsL=sameL&diffL;
                same[t][x]=(byte)BitUtil.pop(sameL);
                diff[t][x]=(byte)BitUtil.pop(diffL);
                hets[t][x]=(byte)BitUtil.pop(hetsL);
            }
        }
        int best=-1;
        for(int x=0; x<sections; x++) {
            ArrayList<Integer> oldCloseHDTaxa=new ArrayList<Integer>();
            for (int t = 0; t < hdTarget.length; t++) {oldCloseHDTaxa.add(t);}
            int window=3;
            while((oldCloseHDTaxa.size()!=1)&&(window<maxWindow)) {
                window++;
                double maxIdentity=0;
                ArrayList<Integer> newCloseHDTaxa=new ArrayList<Integer>();
                for (int t : oldCloseHDTaxa) {
                    double identity=getIdentity(t,x,window);
                    if(Double.isNaN(identity)) {
                        continue;
                    } else if(identity>(maxIdentity+minIdentityDiff)) {
                        newCloseHDTaxa=new ArrayList<Integer>();
                        newCloseHDTaxa.add(t);
                        maxIdentity=identity;
                    } else if(identity>(maxIdentity-minIdentityDiff)) {
                        newCloseHDTaxa.add(t);
                        if(identity>maxIdentity) maxIdentity=identity;
                    }
                }
                if(newCloseHDTaxa.size()>0) oldCloseHDTaxa=newCloseHDTaxa;
            }
            best=oldCloseHDTaxa.get(0);
            int ldAlignTN=hdTarget[best];
//            System.out.printf("For %s section %d window %d best hit %d %d %s \n", ldAlign.getTaxaName(taxa), x, window,
//                    best, ldAlignTN, ldAlign.getTaxaName(ldAlignTN));
//            System.out.print(" "+ldAlign.getTaxaName(ldAlignTN));
            if(x==0) {posBreaks.add(0); hdTaxa.add(best);}
            else if(hdTaxa.get(hdTaxa.size()-1)!=best) {posBreaks.add(ldAlign.getPositionInLocus(x*64)); hdTaxa.add(best);}
        }
//        System.out.println("");
//        System.out.println(ldAlign.getTaxaName(taxa)+"\t"+siteBreaks.size()+"\t"+sections);
        ArrayList[] result={posBreaks, hdTaxa};
        return result;
    }
    
    private double getIdentity(int t, int center, int window) {
        int halfWindow=window-1;
        int sameSum=0, diffSum=0, hetSum=0;
        for(int x=center-halfWindow; x<=center+halfWindow; x++) {
            if((x<0)||(x>=same[0].length)) continue;
            sameSum+=same[t][x];
            diffSum+=diff[t][x];
            hetSum+=hets[t][x];
        }
        int sum=sameSum+diffSum+hetSum;
        int sites=sum-hetSum;
        if(sites<minSites) return Double.NaN;
        double identity=(double)(sameSum+(hetSum/2))/(double)(sum);
        return identity;
    }
    
    private void initHDTargets() {
        hdTarget=new int[hdTargetID.getIdCount()];
        for (int i = 0; i < hdTargetID.getIdCount(); i++) {
            hdTarget[i]=ldAlign.getIdGroup().whichIdNumber(hdTargetID.getIdentifier(i).getFullName());
            System.out.println(hdTargetID.getIdentifier(i).getName()+":"+i+":"+hdTarget[i]);
        }
       
    }
    
    
}
