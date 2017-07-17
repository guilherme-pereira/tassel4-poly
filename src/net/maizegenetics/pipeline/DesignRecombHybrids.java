/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pipeline;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.BitAlignmentHDF5;
import net.maizegenetics.pal.alignment.BitNucleotideAlignment;
import net.maizegenetics.pal.alignment.CombineAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;

/**
 *
 * @author edbuckler
 */
public class DesignRecombHybrids {
 //   private static String root="/Volumes/Thunderbolt SSD/build20120701/06_HapMap/Final/04_BPECFilteredSNPs/";
    private static String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";
    private static String troot="/Users/edbuckler/EdField/RareAllele/B2PDHybrids/LowRecombination/";
    private static String taxaFile="taxaIn2012P.txt";
    
    public static void main(String[] args) {
       // createMergedFiles();
//        BitAlignmentHDF5 nextAlignment=BitAlignmentHDF5.getInstance(root+"comb1to10x.hmp.h5");
//        System.out.println("Open H5 alignment sites:"+nextAlignment.getSiteCount());
//        Alignment fa=AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(nextAlignment, 0.001, 0.999, 4);
//        System.out.println("Filtered alignment"+fa.getSiteCount());
//        Alignment baseA=BitAlignment.getInstance(fa,true);
//        System.out.println("Create memory instance alignment");
//        ExportUtils.writeToHapmap(baseA, false, root+"comb1to10filt.hmp.txt",'\t',null);
        
        Alignment theA=ImportUtils.readFromHapmap(root+"comb1to10filt.hmp.txt",null);
        IdGroup keepIDS=readTaxaList(troot+taxaFile);
        theA=FilterAlignment.getInstance(theA, keepIDS);
        theA.optimizeForTaxa(null);
        OpenBitSet theCentromereMask=makeMask(theA);
        OpenBitSet theArmMask=(OpenBitSet)theCentromereMask.clone();
        theArmMask.not();
        OpenBitSet theAllMask=(OpenBitSet)theCentromereMask.clone();
        theAllMask.union(theArmMask);
        for (int i = 0; i < theA.getSequenceCount(); i++) {
            for (int j = 0; j < i; j++) {
                double centDiv=computeHetBitDistances(theA,i,j,500,theCentromereMask.getBits());
                double armDiv=computeHetBitDistances(theA,i,j,500,theArmMask.getBits());
                double allDiv=computeHetBitDistances(theA,i,j,500,theAllMask.getBits());
                //if(theA.getTaxaName(i).contains("B73")||theA.getTaxaName(i).contains("PI550473")||theA.getTaxaName(i).contains("Ki3")) {
                    System.out.printf("%s %s %g %g %g %n",theA.getTaxaName(i),theA.getTaxaName(j),centDiv,
                            armDiv,allDiv);
               // }
            }
            
        }
        
    }
    
    private static OpenBitSet makeMask(Alignment a) {
        Locus[] l=a.getLoci();
        OpenBitSet obs=new OpenBitSet(a.getSiteCount());
        for(Locus lc: l) {
            int start=lc.getStart();
            int end=lc.getEnd();
            int startSite=Math.abs(a.getSiteOfPhysicalPosition((int)((end-start)*.25), lc));
            int endSite=Math.abs(a.getSiteOfPhysicalPosition((int)((end-start)*.75), lc));
            System.out.printf("%s StartS:%d EndS:%d %n",lc.getChromosomeName(),startSite, endSite);
            for (int i = startSite; i < endSite; i++) {
                obs.fastSet(i); 
            }
            System.out.println("Bits set:"+obs.cardinality());
        }
        return obs;
    } 
    
    private static void createMergedFiles() {
        BitNucleotideAlignment[] eachChr=new BitNucleotideAlignment[10];
        BitAlignmentHDF5 nextAlignment=BitAlignmentHDF5.getInstance(root+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_chr1.hmp.h5");
        IdGroup keepIDS=readTaxaList(root+taxaFile);
        Alignment fa=FilterAlignment.getInstance(nextAlignment, keepIDS);
        Alignment baseA=BitAlignment.getInstance(fa,true);
        baseA.optimizeForTaxa(null);
        baseA.optimizeForSites(null);
        for (int i = 2; i <= 10; i++) {
            System.out.println("Chr"+i);
            nextAlignment=BitAlignmentHDF5.getInstance(root+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_chr"+i+".hmp.h5");
            Alignment a=BitAlignment.getInstance(FilterAlignment.getInstance(nextAlignment, keepIDS),true);
            Alignment[] comb={baseA,a};
            baseA=BitAlignment.getInstance(CombineAlignment.getInstance(comb, true),true);
            baseA.optimizeForTaxa(null);
            baseA.optimizeForSites(null);
        }
        ExportUtils.writeToHDF5(baseA, root+"comb1to10x.hmp.h5");
    }
    
    private static IdGroup readTaxaList(String taxaListFile) {
        List taxa = new ArrayList();
        BufferedReader br = null;
        try {
            br = Utils.getBufferedReader(taxaListFile);
            String inputline = br.readLine();
            Pattern sep = Pattern.compile("\\s+");
            while (inputline != null) {
                inputline = inputline.trim();
                String[] parsedline = sep.split(inputline);
                for (int i = 0; i < parsedline.length; i++) {
                    if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                        taxa.add(parsedline[i]);
                    }
                }
                inputline = br.readLine();
            }
            br.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
        Identifier[] ids = new Identifier[taxa.size()];
        for (int i = 0; i < taxa.size(); i++) {
            ids[i] = new Identifier((String) taxa.get(i));
        }
        return (new SimpleIdGroup(ids));
    }
    
    public static double computeHetBitDistances(Alignment theTBA, int taxon1, int taxon2, 
            int minSitesCompared, long[] mask) {
        
        theTBA = AlignmentUtils.optimizeForTaxa(theTBA);
        long[] iMj = theTBA.getAllelePresenceForAllSites(taxon1, 0).getBits();
        long[] iMn = theTBA.getAllelePresenceForAllSites(taxon1, 1).getBits();
        long[] jMj = theTBA.getAllelePresenceForAllSites(taxon2, 0).getBits();
        long[] jMn = theTBA.getAllelePresenceForAllSites(taxon2, 1).getBits();
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = 0; x < iMj.length; x++) {
            long same = ((iMj[x] & jMj[x]) | (iMn[x] & jMn[x])) & mask[x];
            long diff = ((iMj[x] & jMn[x]) | (iMn[x] & jMj[x])) & mask[x];
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
        }
        double identity = (double) (sameCnt + (hetCnt / 2)) / (double) (sameCnt + diffCnt + hetCnt);
        double dist = 1 - identity;
        int sites = sameCnt + diffCnt - hetCnt;
        if (sites > minSitesCompared) {
            return dist;
        } else {
            return Double.NaN;
        }
        

    }
    
}
