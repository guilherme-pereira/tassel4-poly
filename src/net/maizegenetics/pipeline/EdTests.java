/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pipeline;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.DistanceMatrixUtils;
import net.maizegenetics.pal.ids.*;
import net.maizegenetics.pal.popgen.BitNeighborFinder;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;

/**
 *
 * @author edbuckler
 */
public class EdTests {
 //   String gbsFile="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/allZea20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c10.hmp.txt";
  //  String gbsFile="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c10.hmp.txt";
//    String gbsFile="/Users/edbuckler/SolexaAnal/GBS/build20120110/bpec/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c10.hmp.txt";
    String gbsFile="/Volumes/LaCie/bigprojection/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c10.hmp.txt";
    //String hapFileAGP1="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/chr10.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
//    String hapFileAGP1="/Users/edbuckler/SolexaAnal/HapMapV2/HapMapV2RefgenV2_NAMfounders_teosinte_20110927_chr10.hmp.txt";
//    String hapFileAGP1="/Volumes/LaCie/bigprojection/maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt";
    String hapFileAGP1="/Users/edbuckler/SolexaAnal/bigprojection/maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt";
    String hapFileAGP2="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/chr10.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
    String hapFileAGP3="/Users/edbuckler/SolexaAnal/HapMapV2/jermpipe/SNPS201010/fusion/chr10.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
    
    String gbsHapMergeFile="/Volumes/LaCie/bigprojection/MergeGBSHap.c10.hmp.txt";
 //   String gbsHapMergeImpFile="/Volumes/LaCie/bigprojection/MergeGBSHapImp.c10.hmp.txt";
    String gbsHapMergeImpFile="/Users/edbuckler/SolexaAnal/bigprojection/MergeGBSHapImp.c10.hmp.txt";
    BitAlignment gbsMap=null;
    BitAlignment hapMap=null;
    BitAlignment mergeMap=null;
    

    public EdTests() {
//        String header="/Users/edbuckler/SolexaAnal/bigprojection/";
        String header="/Volumes/LaCie/bigprojection/";
//        this.createProjAlignment(header+"maizeHapMapV2_B73RefGenV2_201203028_chr4.hmp.txt", 
//                header+"Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c4.hmp.txt", header+"chr4.projA.txt", 4);
        int c=6;
        this.createProjAlignment(header+"maizeHapMapV2_B73RefGenV2_201203028_chr"+c+".hmp.txt", 
                header+"Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c"+c+".hmp.txt", header+"chr"+c+".imp95_1024.projA.txt", 4);
//        this.createProjAlignment(header+"maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt", 
//                header+"282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c10.hmp.txt", header+"282chr10.projA.txt", 4);
        
//        this.exportRegion(header+"282chr10.projA.txt", header+"282chr10_94M_95M.hmp.txt", 94000000,95000000);
        
//
        int s=81, e=83;
        this.exportRegion(header+"chr"+c+".imp95_1024.projA.txt", header+"chr"+c+".imp95_"+s+"M_"+e+"M.hmp.txt",
                "/Volumes/LaCie/bigprojection/amestaxa.txt",
                s*1000000,e*1000000);
        s=83;
        e=85;
        this.exportRegion(header+"chr"+c+".imp95_1024.projA.txt", header+"chr"+c+".imp95_"+s+"M_"+e+"M.hmp.txt",
                "/Volumes/LaCie/bigprojection/amestaxa.txt",
                s*1000000,e*1000000);
//        SBitAlignment hapMapOld=(SBitAlignment)ImportUtils.readFromHapmap(hapFileAGP3, (ProgressListener)null);
//        printTaxaNames(hapMapOld);
        System.exit(0);
        
 //       convertFilesToFast(false, true);
 //       System.exit(0);
//        gbsMap=(TBitAlignment)readGZOfSBit(gbsFile, false);
//        hapMap=(SBitAlignment)readGZOfSBit(hapFileAGP1, true);
//        hapMap=(SBitAlignment)fixHapMapNames(hapMap);  //adds tags so that HapMapNames are recognizable
//        Alignment mna=combineAlignments(hapMap, gbsMap);
//        ExportUtils.writeToHapmap(mna, true, gbsHapMergeFile, '\t', null);
        
//        SBitAlignment temp=(SBitAlignment)ImportUtils.readFromHapmap(gbsHapMergeFile, (ProgressListener)null);
//        mergeMap=TBitAlignment.getInstance(temp);
//        Alignment mna=imputeHapMapTaxaWithNearIdenticals(mergeMap);
//        ExportUtils.writeToHapmap(mna, true, gbsHapMergeImpFile, '\t', null);
//        System.exit(0);
//        convertFilesToFastTbit(gbsHapMergeImpFile);
        gbsMap=(BitAlignment)readGZOfSBit(gbsHapMergeImpFile, false);
        hapMap=(BitAlignment)readGZOfSBit(hapFileAGP1, true);
        hapMap=(BitAlignment)fixHapMapNames(hapMap); 
//        compareIdentity("B73", hapMap, "B73", gbsMap, true);
//        compareIdentity("B73", hapMap, "B97", gbsMap, true);
//        compareIdentity("B97", hapMap, "B97", gbsMap, true);
//        compareIdentity("B97", hapMap, "B73", gbsMap, true);
//        compareIdentity("CML277", hapMap, "CML277", gbsMap, true);
//        compareIdentity("CML277", hapMap, "B73", gbsMap, true);
//        compareIdentity("TIL08", hapMap, "TIL08", gbsMap, true);
//        System.exit(0);
        System.out.println("GBS Map Taxa");
        printTaxaNames(gbsMap);
        System.out.println("HapMap Taxa");
        printTaxaNames(hapMap);
//        for(int i=0; i<hapMap.getSequenceCount(); i++) {
//            System.out.print("\""+hapMap.getTaxaName(i)+"\",");
//        }
//        System.out.println("");
       String[] taxa={"B73","B97","CML103","CML228","CML247","CML277","CML322","CML333","CML52","CML69","HP301",
           "IL14H","Ki11","Ki3","KY21","M162W","M37W","MO17","MO18W","MS71","NC350","NC358","OH43","OH7B","P39",
           "TIL01","TIL02","TIL03","TIL04","TIL05","TIL06","TIL07","TIL08","TIL09","TIL10","TIL11","TIL12",
           "TIL15","TIL16","TIL17","TIL25","TIL14","TX303","TZI8","W22","W64A"};
       IdGroup testHD=new SimpleIdGroup(taxa);
       //I am using this testHD as the taxa names from hapMap do not match gbsMap.  
       //We will probably need to load a TreeMap to redirect between the two.
       BitNeighborFinder bnf=new BitNeighborFinder(hapMap.getIdGroup(), gbsMap, hapMap);
       Alignment pa=bnf.getPa();
       compareIdentity("B73", hapMap, "Z001E0128", pa, false);
       compareIdentity("B73", hapMap, "Z001E0101", pa, false);
       compareIdentity("B97", hapMap, "Z001E0101", pa, false);
       compareIdentity("B73", hapMap, "M0236", pa, false);
       compareIdentity("B73", hapMap, "W22", pa, false);
        System.out.println("Within HapMap");
        compareIdentity("B73", hapMap, "B73", hapMap, false);
       compareIdentity("B73", hapMap, "W22", hapMap, false);
       compareIdentity("B73", hapMap, "MO17", hapMap, false);
       compareIdentity("B73", hapMap, "B97", hapMap, false);
       System.out.println("Within PA");
       compareIdentity("B73(PI550473)", pa, "Z001E0128", pa, false);
       compareIdentity("B73(PI550473)", pa, "Z001E0101", pa, false);
       compareIdentity("B97", pa, "Z001E0101", pa, false);
       compareIdentity("B73(PI550473)", pa, "M0236", pa, false);
       compareIdentity("B73(PI550473)", pa, "W22", pa, false);
       
        
        
//        compareSitesInFiles();
 //       hapMap2=new MutableTBitNucleotideAlignment(hapMap, hapMap.getMaxNumAlleles(), hapMap.retainsRareAlleles());
//        convertHapMapSites();
//        hapMap=hapMap2;
//        compareSitesInFiles();
    }
    
    public void createProjAlignment(String hFile, String gFile, String pOutFile, int chr) {
        Alignment gbsMap = ImportUtils.readFromHapmap(gFile, false, null);
        //BitAlignment gbsMap=TBitAlignment.getInstance(ImportUtils.readFromHapmap(gFile, (ProgressListener)null));
        System.out.println("GBS Map Read");
//        SBitAlignment hapMap=(SBitAlignment)readGZOfSBit(hapFileAGP1, true);
        Alignment hapMap=ImportUtils.readFromHapmap(hFile, true, (ProgressListener)null);
        System.out.println("HapMap Read");
        hapMap=fixHapMapNames(hapMap);  //adds tags so that HapMapNames are recognizable
        System.out.println("HapMap Names Fixed");
        MutableNucleotideAlignment mna=combineAlignments(hapMap, gbsMap);
        System.out.println("HapMap and GBS combined");
        mna.clean();
        Alignment mergeMap=BitAlignment.getInstance(mna, false);
        System.out.println("Merge convertd to TBit");
        mna=imputeHapMapTaxaWithNearIdenticals(mergeMap);
        mna.clean();
        Alignment gbsImpMap=BitAlignment.getInstance(mna, false);
        System.out.println("Imputed HapMapTaxa on MergeMap");
        BitNeighborFinder bnf=new BitNeighborFinder(hapMap.getIdGroup(), gbsImpMap, hapMap);
        ProjectionAlignment pa=bnf.getPa();
        pa.reportPAComposition();
        writeAlignmentToSerialGZ(pa, pOutFile);
    }
    
    public void exportRegion(String paFile, String outFile, String taxaListFile, int start, int end) {
        ProjectionAlignment pa=readGZOfPA(paFile);
//        pa.reportPAComposition();
        System.out.println("Done reading Projection Alignment");
        int startSite=Math.abs(pa.getSiteOfPhysicalPosition(start, null));
        int endSite=Math.abs(pa.getSiteOfPhysicalPosition(end, null));
        Alignment fa=FilterAlignment.getInstance(pa, startSite, endSite);
        System.out.printf("Filtering Complete startSite:%d endSite:%d %n", startSite, endSite);
        if(taxaListFile!=null) {
            IdGroup subIDS=getIDToKeepFromFile( taxaListFile);
            fa=FilterAlignment.getInstance(fa, subIDS,false);
            System.out.printf("Filtering Complete taxa:%d %n", fa.getSequenceCount());
        }
        
        ExportUtils.writeToHapmap(fa, false, outFile, '\t', null);
        System.out.println("Export complete");
    }
    
    public static Alignment fixHapMapNames(Alignment hapMap) {
        for (int i = 0; i < hapMap.getIdGroup().getIdCount(); i++) {
            Identifier id=hapMap.getIdGroup().getIdentifier(i);
            String s=id.getFullName();
            s=s+":HMP";
            s=s.replaceFirst("KI3", "Ki3");
            s=s.replaceFirst("KI11", "Ki11");
            s=s.replaceFirst("CML312SR", "CML312");
            s=s.replaceFirst("W22", "W22q");
            s=s.replaceFirst("TIL04-TIP454", "TIL04");
            hapMap.getIdGroup().setIdentifier(i, new Identifier(s));
        }
        return hapMap;
    }
    
    private int[] indicesOfHapMapTaxa(IdGroup hapIDs) {
        ArrayList<Integer> hid=new ArrayList<Integer>();
        for (int i = 0; i < hapIDs.getIdCount(); i++) {
            if(hapIDs.getIdentifier(i).getFullName().contains("HMP")) hid.add(i);
        }
        int[] ahid=new int[hid.size()];
        for (int i = 0; i < ahid.length; i++) {
            ahid[i]=hid.get(i);
        }
        return ahid;
    }
    
    private MutableNucleotideAlignment imputeHapMapTaxaWithNearIdenticals(Alignment mergeAlign) {
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(mergeAlign);
        int[] hapMapTaxaInd=indicesOfHapMapTaxa(mergeAlign.getIdGroup());
        for (int t = 0; t < hapMapTaxaInd.length; t++) {
            int ht=hapMapTaxaInd[t];
            TreeMap<Double,Integer> closestTaxon=new TreeMap<Double,Integer>();
            for (int gt = 0; gt < mergeAlign.getSequenceCount(); gt++) {
                if(ht==gt) continue;
                double dist=DistanceMatrixUtils.getIBSDistance(mergeAlign.getAllelePresenceForAllSites(gt, 0), 
                        mergeAlign.getAllelePresenceForAllSites(gt, 1), mergeAlign.getAllelePresenceForAllSites(ht, 0),
                        mergeAlign.getAllelePresenceForAllSites(ht, 1));
                if(dist<0.05) closestTaxon.put(dist, gt);
            }
            System.out.println(mergeAlign.getTaxaName(ht));
            System.out.println(closestTaxon.toString());
            if(closestTaxon.size()<1) continue;
            ArrayList<Integer> taxaList=new ArrayList(closestTaxon.values());
            for (int c : taxaList) {
              for (int s = 0; s < mna.getSiteCount(); s++) {
                  if(mna.getBase(ht, s)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                      if(!AlignmentUtils.isHeterozygous(mna.getBase(c,s))) mna.setBase(ht, s, mna.getBase(c,s));
                  }
              }
            }
        }
        return mna;
    }
    
    public static MutableNucleotideAlignment combineAlignments(Alignment hapMap, Alignment gbsAlign) {
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(gbsAlign,
                hapMap.getSequenceCount()+gbsAlign.getSequenceCount(), gbsAlign.getSiteCount());
        System.out.println("MNA Created");
        int orgTaxaNum=mna.getSequenceCount();
        for (int i = 0; i < hapMap.getSequenceCount(); i++) {
            Identifier id=hapMap.getIdGroup().getIdentifier(i);
            mna.addTaxon(id);  
        }
        mna.clean();
        System.out.println("MNA taxa added and cleaned");
        int mappingSites=0, sameStateSites=0;
        for (int i = 0; i < hapMap.getSiteCount(); i++) {
            int pos=hapMap.getPositionInLocus(i);
            int site=mna.getSiteOfPhysicalPosition(pos, null);
            if(site<0) continue;
            mappingSites++;
            if(hapMap.getMajorAllele(i)!=mna.getMajorAllele(site)) continue;
            if(hapMap.getMinorAllele(i)!=mna.getMinorAllele(site)) continue;
            sameStateSites++;
            for (int j = 0; j < hapMap.getSequenceCount(); j++) {
                byte hp=hapMap.getBase(j, i);
                if (AlignmentUtils.isHeterozygous(hp)) continue;
                mna.setBase(orgTaxaNum+j, site, hapMap.getBase(j, i));
            }
        }
        System.out.printf("Mapping Sites=%d SameStates=%d %n", mappingSites, sameStateSites);
        System.out.println("MNA taxa bases added");
        int[] siteComp=new int[mna.getSiteCount()];
        int[] siteErrors=new int[mna.getSiteCount()];
        int[] taxaComp=new int[mna.getSiteCount()];
        int[] taxaErrors=new int[mna.getSiteCount()];
        for (int j = 0; j < hapMap.getSequenceCount(); j++) {          
            int tgbs=gbsAlign.getIdGroup().whichIdNumber(hapMap.getTaxaName(j));
            if(tgbs<0) continue;
            System.out.println(mna.getTaxaName(tgbs)+":"+mna.getTaxaName(j+orgTaxaNum));
            for (int s = 0; s < mna.getSiteCount(); s++) {
                byte hb=mna.getBase(j+orgTaxaNum, s);
                byte gb=mna.getBase(tgbs, s);
                if(((hb >>> 4) & 0xf)!=(hb & 0xf)) continue;  //heterozgyous test
                if(((gb >>> 4) & 0xf)!=(gb & 0xf)) continue;
                if(hb==Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                    mna.setBase(j, s, gb);
                }  else if(gb!=Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    siteComp[s]++;
                    taxaComp[j]++;
                    if(hb!=gb) {
                        siteErrors[s]++;
                        taxaErrors[j]++;
  //                      System.out.println(NucleotideAlignmentConstants.getNucleotideIUPAC(hb)+":"+NucleotideAlignmentConstants.getNucleotideIUPAC(gb));
                    }
                }
            }
        }
        System.out.println("Merged when possible");
        for (int j = 0; j < hapMap.getSequenceCount(); j++) {
            System.out.printf("s%s %d %d %n",hapMap.getTaxaName(j),taxaComp[j], taxaErrors[j]);
        }
//        for (int i = 0; i < siteComp.length; i++) {
//            System.out.printf("s%d %d %d %n",i,siteComp[i], siteErrors[i]);
//        }
        return mna;
    }
    
    
    
    
    public void printTaxaNames(Alignment a) {
        IdGroup idg=a.getIdGroup();
        for (int i = 0; i < idg.getIdCount(); i++) {
            System.out.printf("%d \t %s \t %s %n",i,idg.getIdentifier(i).getName(),idg.getIdentifier(i).getFullName());     
        }
    }
    
    private void compareIdentity(String taxon, Alignment hapMap, String taxon2, Alignment projAlign, boolean withSiteLookup) {
        int t1=hapMap.getIdGroup().whichIdNumber(taxon);
        if(hapMap instanceof ProjectionAlignment) {
            System.out.println(taxon+">"+((ProjectionAlignment)hapMap).getCompositionOfTaxon(t1));
        }
        int t2=projAlign.getIdGroup().whichIdNumber(taxon2);
        if(projAlign instanceof ProjectionAlignment) {
            System.out.println(taxon2+">"+((ProjectionAlignment)projAlign).getCompositionOfTaxon(t2));
        }
        int same=0, diff=0;
        for (int i = 0; i < hapMap.getSiteCount(); i++) {
            if(hapMap.isHeterozygous(t1, i)) continue;
            if((hapMap.getMajorAllele(i)==4)||(hapMap.getMajorAllele(i)==5)) continue;
            if((hapMap.getMinorAllele(i)==4)||(hapMap.getMinorAllele(i)==5)) continue;
            byte b1=hapMap.getBase(t1, i);
            byte b2=Alignment.UNKNOWN_DIPLOID_ALLELE;
            int site2=i;
            if(withSiteLookup) {
                int pos=hapMap.getPositionInLocus(i);
                site2=projAlign.getSiteOfPhysicalPosition(pos, null);
                if(site2>-1) b2=projAlign.getBase(t2, site2);
                
            } else {
                b2=projAlign.getBase(t2, site2);
            }
//            if(projAlign.isHeterozygous(t2, site2)) continue;
            if(b1==Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
            if(b2==Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
            if(b1==b2) {same++;}
            else {diff++;
//                System.out.println(hapMap.getBaseAsString(t1, i)+":"+projAlign.getBaseAsString(t2, site2));
            }
//            if(i%100000==0) System.out.printf("%d %s  %s  %d %d %g %n",i, taxon, taxon2, same, diff, ((double)same/(double)(same+diff)));
        }
        System.out.printf("%s  %s  %d %d %g %n",taxon, taxon2, same, diff, ((double)same/(double)(same+diff)));
    }
  
    
    private void convertFilesToFast(boolean gbsSbit, boolean hapSbit) {
        Alignment gbs=ImportUtils.readFromHapmap(gbsFile, true, (ProgressListener)null);
        if(gbsSbit) {writeAlignmentToSerialGZ(gbs, gbsFile);}
        else{ 
            Alignment gbsOut=BitAlignment.getInstance(gbs, false);
            writeAlignmentToSerialGZ(gbsOut, gbsFile);
        }
        Alignment hapmapIn=ImportUtils.readFromHapmap(hapFileAGP1, true, (ProgressListener)null);
        if(hapSbit) {writeAlignmentToSerialGZ(hapmapIn, hapFileAGP1);}
        else {
            Alignment hapmapOut=BitAlignment.getInstance(hapmapIn, false);
            writeAlignmentToSerialGZ(hapmapOut, hapFileAGP1);
        }
    }
    
    private void convertFilesToFastSbit(String flatFile) {
        Alignment align=ImportUtils.readFromHapmap(flatFile, true, (ProgressListener)null);
        writeAlignmentToSerialGZ(align, flatFile);
    }
    
    private void convertFilesToFastTbit(String flatFile) {
        Alignment align=ImportUtils.readFromHapmap(flatFile, true, (ProgressListener)null);
        Alignment alignT=BitAlignment.getInstance(align, false);
        writeAlignmentToSerialGZ(alignT, flatFile);
    }
    
    private void compareSitesInFiles() {
        int count=0;
        Locus lx=hapMap.getLocus(0);
        for (int i = 0; i < gbsMap.getSiteCount(); i++) {
            int position=gbsMap.getPositionInLocus(i);
            int hapSite=hapMap.getSiteOfPhysicalPosition(position, lx);
            if(hapSite>0) System.out.printf("%d %d %d %n",count++, position, hapSite);
        }
    }
    
     
    public static void main(String[] args) {
        EdTests et=new EdTests();
        
    }
    
    static ProjectionAlignment readGZOfPA(String inFile) {
        ProjectionAlignment pa=null;
        long time=System.currentTimeMillis();
        try {
            File theFile = new File(Utils.addSuffixIfNeeded(inFile, ".gz"));
            System.out.println("Reading PA:"+theFile);
            FileInputStream fis = new FileInputStream(theFile);
            GZIPInputStream gs = new GZIPInputStream(fis);
            ObjectInputStream ois = new ObjectInputStream(gs);
            pa=(ProjectionAlignment)ois.readObject();
            System.out.println(pa.getSiteCount());
            System.out.println(pa.getSequenceCount());
            ois.close();
            fis.close();
            
        } catch (Exception ee) {
            //sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }
        System.out.println("Time:"+(System.currentTimeMillis()-time));
        return pa;
    }
    
    
    static Alignment readGZOfSBit(String inFile, boolean isSBit) {
        BitAlignment sba=null;
        BitAlignment tba=null;
        long time=System.currentTimeMillis();
        try {
            File theFile = new File(Utils.addSuffixIfNeeded(inFile, ".gz"));
            System.out.println("Reading:"+theFile);
            FileInputStream fis = new FileInputStream(theFile);
            GZIPInputStream gs = new GZIPInputStream(fis);
            ObjectInputStream ois = new ObjectInputStream(gs);
            if(isSBit) {sba=(BitAlignment)ois.readObject();
                System.out.println(sba.getSiteCount());
                System.out.println(sba.getSequenceCount());}
            else {tba=(BitAlignment)ois.readObject();
                System.out.println(tba.getSiteCount());
                System.out.println(tba.getSequenceCount());}
            ois.close();
            fis.close();
            
        } catch (Exception ee) {
            //sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }
        System.out.println("Time:"+(System.currentTimeMillis()-time));
        if(isSBit) return sba;
        return tba;
    }
    
    static void writeAlignmentToSerialGZ(Alignment sba, String outFile) {
        //String inFile="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/allZea20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c10.hmp.txt";
        long time=System.currentTimeMillis();
        try {

            File theFile = new File(Utils.addSuffixIfNeeded(outFile, ".gz"));
            FileOutputStream fos = new FileOutputStream(theFile);
            GZIPOutputStream gz = new GZIPOutputStream(fos);
            ObjectOutputStream oos = new ObjectOutputStream(gz);
            oos.writeObject(sba);
            oos.flush();
            oos.close();
            fos.close();
            System.out.println("Wrote:"+theFile.toString());
        } catch (Exception ee) {
            //sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }
        System.out.println("Time:"+(System.currentTimeMillis()-time));
    }
    
    static IdGroup getIDToKeepFromFile(String taxaListFile) {
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
        } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        Identifier[] ids = new Identifier[taxa.size()];
        for (int i = 0; i < taxa.size(); i++) {
            ids[i] = new Identifier((String) taxa.get(i));
        }
        return (new SimpleIdGroup(ids));
    }
}
