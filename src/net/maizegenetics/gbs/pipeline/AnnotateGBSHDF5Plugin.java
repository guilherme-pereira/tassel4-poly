/*
 * BiParentalErrorCorrectionPlugin
 */
package net.maizegenetics.gbs.pipeline;


import cern.colt.list.IntArrayList;
import ch.systemsx.cisd.hdf5.HDF5DataClass;
import ch.systemsx.cisd.hdf5.HDF5DataSetInformation;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.SiteMappingInfo;
import net.maizegenetics.gwas.imputation.NucleotideImputationUtils;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.BitAlignmentHDF5;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.IBSErrorByTaxon;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.pal.statistics.FisherExact;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import org.apache.log4j.Logger;

/**
 * Tools for characterizing all the SNPs with MAF, LD patterns, mapping distribution, etc.
 *
 * Error rates are bounded away from zero, but adding 0.5 error to all error
 * rates that that were observed to be zero.
 *
 * @author edbuckler
 */
public class AnnotateGBSHDF5Plugin extends AbstractPlugin {
    private int minSite=500;
    private double IBSthreshold=0.01;
    private float[] errorRate, errorsMinorRate;
    //should we focus on minor allele agreement??
    private int[] errorCnt, mjCorrCnt, mnCorrCnt;
    private double maxErrorRate = 0.05;
    private boolean myRemoveUntestedError = true;
    private String outHapMap = null;
    private static ArgsEngine engine = new ArgsEngine();
    int start = 1, end = 1;
    private String infile, outfileAnnoHDF5;
    private static final Logger myLogger = Logger.getLogger(AnnotateGBSHDF5Plugin.class);
    int homoweight = 1;  //weight of homozygous genotypes to heterozgyous genotypes
    //1 assumes low coverage and each homozygote is only counted once, while hets get a count for each allele
    //2 assumes better coverage, and homozygotes are counted twice.
    private String snpLogFileName;
    private SNPLogging snpLogging = null;
    private byte[] minorGenotype;

    public AnnotateGBSHDF5Plugin() {
        super(null, false);
    }

    public AnnotateGBSHDF5Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }
  
    public void reportPercentilesOfErrorRates() {
        //TODO make it saveable to file or stdout
        float[] le = Arrays.copyOf(errorRate, errorRate.length);
        Arrays.sort(le);
        System.out.println("Percentile\tErrorRate");
        for (int i = 0; i < le.length; i += (le.length / 20)) {
            System.out.printf("%.2g\t%.3g%n", ((double) i / (double) le.length), le[i]);
        }
    }
    
   public static void reportPercentilesOfErrorRates(double[] arr, int intervals) {
        //TODO make it saveable to file or stdout
        double[] le = Arrays.copyOf(arr, arr.length);
        Arrays.sort(le);
        System.out.println("Percentile\tErrorRate");
        for (int i = 0; i < le.length; i += (le.length / intervals)) {
            System.out.printf("%.2g\t%.3g%n", ((double) i / (double) le.length), le[i]);
        }
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-oE", "--outErrorFile", true);
        engine.add("-oB", "--outBinDistFile", true);
        engine.add("-sC", "--startChrom", true);
        engine.add("-eC", "--endChrom", true);
        engine.add("-mxE", "--maxError", true);
        engine.add("-kpUT", "--keepUntested", false);
        engine.add("-snpLog", "", true);
        engine.parse(args);
        if (engine.getBoolean("-sC")) {
            start = Integer.parseInt(engine.getString("-sC"));
        }
        if (engine.getBoolean("-eC")) {
            end = Integer.parseInt(engine.getString("-eC"));
        }
        infile = engine.getString("-hmp");
        outfileAnnoHDF5 = engine.getString("-o");
        if (engine.getBoolean("-snpLog")) {
            snpLogFileName = engine.getString("-snpLog");
        }
        snpLogging = new SNPLogging(snpLogFileName, this.getClass());
        performFunction(null);
    }

    public void setOutHapMap(String outHapMap) {
        this.outHapMap = outHapMap;
    }

    public void setMaxErrorRate(double maxErrorRate) {
        this.maxErrorRate = maxErrorRate;
    }


    private void printUsage() {
        myLogger.info(
                "\n\n\nAvailable options for the BiParentalErrorCorrectionPlugin are as follows:\n"
                + "-hmp   Input HapMap file\n"
                + "-o     Output HapMap file\n"
                + "-sC    Start chromosome\n"
                + "-eC    End chromosome\n"
                + "-mxE    Maximum error\n"
                + "-kpUT   Keep untested SNPs for error (default remove)\n"
                + "-snpLog  SNPs Removed Log file name\n\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        while (start <= end) {
            int chr = start;
            ArrayList<Datum> dList = new ArrayList<Datum>();
            String currFile = infile.replace("+", "" + chr);
            System.out.println("Reading: " + currFile);
            BitAlignmentHDF5 align;
            Alignment a;
            try {
             //   align = BitAlignmentHDF5.getInstance(currFile);
                a=ImportUtils.readFromHapmap(infile, null);

            } catch (Exception e) {
                myLogger.info("Could not read input hapmap file for chr" + chr + ":\n\t" + currFile + "\n\tSkipping...");
                continue;
            }
            System.out.println("Finished Reading: " + currFile);
            String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";
//            String outfile=root+"annoNoHets.hmp.h5";
//            annotateMAF(a,  outfileAnnoHDF5);
            annotateIBSError(a,  outfileAnnoHDF5);
//            annotateLD(a, outfileAnnoHDF5,128,4,0.05f, 5000);
            start++;
        }
        snpLogging.close();
        return null;
    }
    
    private void annotateMAF(Alignment ah5, String hdf5File) {
        float[] maf=new float[ah5.getSiteCount()];
        float[] hets=new float[ah5.getSiteCount()];
        float[] scov=new float[ah5.getSiteCount()];
        float taxa=(float)ah5.getSequenceCount();
        for (int i = 0; i < maf.length; i++) {
            maf[i]=(float)ah5.getMinorAlleleFrequency(i);
            int siteCovCnt=ah5.getTotalGametesNotMissing(i)/2;
            scov[i]=(float)siteCovCnt/taxa;
            hets[i]=(float)ah5.getHeterozygousCount(i)/siteCovCnt;  //terry's proportion hets is off all sites
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        myLogger.info("Annotating HDF5 file: " + hdf5File);
//        config.overwrite();
//        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        if(!h5w.exists(HapMapHDF5Constants.SITE_DESC)) h5w.createGroup(HapMapHDF5Constants.SITE_DESC);
        if(!h5w.exists(HapMapHDF5Constants.MAF)) h5w.createFloatArray(HapMapHDF5Constants.MAF, maf.length);
        h5w.writeFloatArray(HapMapHDF5Constants.MAF, maf);
        if(!h5w.exists(HapMapHDF5Constants.SITEHET)) h5w.createFloatArray(HapMapHDF5Constants.SITEHET, hets.length);
        h5w.writeFloatArray(HapMapHDF5Constants.SITEHET, hets);
        if(!h5w.exists(HapMapHDF5Constants.SITECOV)) h5w.createFloatArray(HapMapHDF5Constants.SITECOV, scov.length);
        h5w.writeFloatArray(HapMapHDF5Constants.SITECOV, scov);
    }
    
    
    /**
     * Methods to evaluate the LD of a site.  Currently it evaluates everything by r2 and p-value.  It may be
     * reasonable to rescale the p-value by comparison to distant locations in the genome
     * @param a
     * @param hdf5File
     * @param minPhysDist
     * @param minMinorCnt
     * @param minR2
     * @param numberOfTests 
     */
    private void annotateLD(Alignment a, String hdf5File, int minPhysDist, int minMinorCnt, 
            float minR2, int numberOfTests) {
        FisherExact myFisherExact=new FisherExact(a.getSequenceCount() + 10);
        a.optimizeForSites(null);
        int minCnt=40;
        int sites=a.getSiteCount();
        int binOf10=8;
        float[] maxR2=new float[sites];  Arrays.fill(maxR2, Float.NaN);
        float[] minPOfMaxR2=new float[sites];  Arrays.fill(minPOfMaxR2, Float.NaN);
        int[] minSigDistance=new int[sites]; Arrays.fill(minSigDistance, Integer.MAX_VALUE);
        float[] propSigTests=new float[sites]; Arrays.fill(propSigTests, Float.NaN);
        double maxSum=0, sumSites=0;
        for (int i = 0; i < a.getSiteCount(); i++) {
            int[][] contig = new int[2][2];
            Locus myLocus=a.getLocus(i);
            BitSet rMj = a.getAllelePresenceForAllTaxa(i, 0);
            BitSet rMn = a.getAllelePresenceForAllTaxa(i, 1); 
            if(rMn.cardinality()<minMinorCnt) continue;
            TreeMap<Double,SiteMappingInfo> bestLDSites=new TreeMap<Double,SiteMappingInfo>(Collections.reverseOrder());
            double minLD=0;
            int attemptTests=0, completedTests=0, sigTests=0;
            int leftSite=i, rightSite=i, j=-1;
            int position=a.getPositionInLocus(i), dist=-1, minSigDist=Integer.MAX_VALUE;
            while((completedTests<numberOfTests)&&((leftSite>0)||(rightSite+1<sites))) {
                int rightDistance=(rightSite+1<sites)?a.getPositionInLocus(rightSite+1)-position:Integer.MAX_VALUE;
                int leftDistance=(leftSite>0)?position-a.getPositionInLocus(leftSite-1):Integer.MAX_VALUE;
                if(rightDistance<leftDistance) {rightSite++; j=rightSite; dist=rightDistance;}
                else {leftSite--; j=leftSite; dist=leftDistance;}
                if(dist<minPhysDist) continue;
                attemptTests++;
               // System.out.printf("Dist: %d %d bin: %d %n",dist, Integer.highestOneBit(dist), bin);
                BitSet cMj = a.getAllelePresenceForAllTaxa(j, 0);  //major alleles
                BitSet cMn = a.getAllelePresenceForAllTaxa(j, 1);  //minor alleles
                int n = 0;
                n += contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);
                n += contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
                if(contig[1][0]+contig[1][1]<minMinorCnt) continue;
                n += contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
                if(contig[0][1]+contig[1][1]<minMinorCnt) continue;
                n += contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
                if(n<minCnt) continue;
                double rValue = LinkageDisequilibrium.calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minCnt);
                if (Double.isNaN(rValue)) continue;
                completedTests++;
                if(rValue<minR2) continue;
                double pValue=myFisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
                if(pValue>(0.05/numberOfTests)) continue;
                if(dist<minSigDist) minSigDist=dist;
                sigTests++;
                if(rValue>minLD) {
                    float[] result={j, (float)rValue, (float)pValue, dist};
                    SiteMappingInfo smi=new SiteMappingInfo(Integer.parseInt(a.getLocusName(j)),
                            (byte)1,a.getPositionInLocus(j),(float)rValue, (float)pValue,j);
                    bestLDSites.put(rValue, smi);
                }
                if(bestLDSites.size()>20) {
                    bestLDSites.remove(bestLDSites.lastKey());
                    minLD=bestLDSites.lastKey();
                }
            }
            
            if(bestLDSites.size()>0) {
                SiteMappingInfo smi=(SiteMappingInfo)bestLDSites.firstEntry().getValue();
                maxR2[i]=smi.r2;
                minPOfMaxR2[i]=smi.mapP;
                maxSum+=bestLDSites.firstEntry().getKey();
                sumSites+=bestLDSites.size();
            }
            minSigDistance[i]=minSigDist;
            propSigTests[i]=(float)sigTests/(float)completedTests;
            System.out.printf("s:%d Attempted:%d Completed:%d SigTests:%d MinDist:%d TreeSize: %d %n",i,attemptTests,completedTests,
                    sigTests, minSigDist, bestLDSites.size());
//            System.out.printf("s:%d %s %s %s %n",i,Arrays.toString(testsLD[i]),Arrays.toString(maxR2[i]),Arrays.toString(minP[i]));
            System.out.printf("s:%d %s %n",i,bestLDSites.toString());
            
            if(i%1000==0) System.out.println("Avg MaxR2:"+(maxSum/(double)i)+"  avgSites:"+(sumSites/(double)i));
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        myLogger.info("Annotating HDF5 file with LD: " + hdf5File);
//        config.overwrite();
//        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        if(!h5w.exists(HapMapHDF5Constants.LD_DESC)) h5w.createGroup(HapMapHDF5Constants.LD_DESC);
        if(!h5w.exists(HapMapHDF5Constants.LDR2_DESC)) h5w.createFloatArray(HapMapHDF5Constants.LDR2_DESC, maxR2.length);
        h5w.writeFloatArray(HapMapHDF5Constants.LDR2_DESC, maxR2);
        if(!h5w.exists(HapMapHDF5Constants.LDP_DESC)) h5w.createFloatArray(HapMapHDF5Constants.LDP_DESC, minPOfMaxR2.length);
        h5w.writeFloatArray(HapMapHDF5Constants.LDP_DESC, minPOfMaxR2);
        if(!h5w.exists(HapMapHDF5Constants.LDPropLD_DESC)) h5w.createFloatArray(HapMapHDF5Constants.LDPropLD_DESC, propSigTests.length);
        h5w.writeFloatArray(HapMapHDF5Constants.LDPropLD_DESC, propSigTests);
        if(!h5w.exists(HapMapHDF5Constants.LDMinDist_DESC)) h5w.createIntArray(HapMapHDF5Constants.LDMinDist_DESC, minSigDistance.length);
        h5w.writeIntArray(HapMapHDF5Constants.LDMinDist_DESC, minSigDistance);
        
    }
    
   
    

    private void annotateIBSError(Alignment a, String hdf5File) {
        a.optimizeForTaxa(null);
        int sites=a.getSiteCount();
        errorCnt=new int[sites];
        mjCorrCnt=new int[sites];
        mnCorrCnt=new int[sites];
        for (int bt = 0; bt < a.getSequenceCount(); bt++) {
            IBSErrorByTaxon iebt=new IBSErrorByTaxon(bt,a,75, 2000, 400,0.01);
            mjCorrCnt=addTwoVector(mjCorrCnt,iebt.getMajorCorrectCounts());
            mnCorrCnt=addTwoVector(mnCorrCnt,iebt.getMinorCorrectCounts());
            errorCnt=addTwoVector(errorCnt,iebt.getErrorCounts());
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        myLogger.info("Annotating HDF5 file: " + hdf5File);
//        config.overwrite();
//        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        if(!h5w.exists(HapMapHDF5Constants.SITE_DESC)) h5w.createGroup(HapMapHDF5Constants.SITE_DESC);
        if(!h5w.exists(HapMapHDF5Constants.IBSMAJORCORR_DESC)) h5w.createIntArray(HapMapHDF5Constants.IBSMAJORCORR_DESC, mjCorrCnt.length);
        h5w.writeIntArray(HapMapHDF5Constants.IBSMAJORCORR_DESC, mjCorrCnt);
        if(!h5w.exists(HapMapHDF5Constants.IBSMINORCORR_DESC)) h5w.createIntArray(HapMapHDF5Constants.IBSMINORCORR_DESC, mnCorrCnt.length);
        h5w.writeIntArray(HapMapHDF5Constants.IBSMINORCORR_DESC, mnCorrCnt);
        if(!h5w.exists(HapMapHDF5Constants.IBSERROR_DESC)) h5w.createIntArray(HapMapHDF5Constants.IBSERROR_DESC, errorCnt.length);
        h5w.writeIntArray(HapMapHDF5Constants.IBSERROR_DESC, errorCnt);
        float[] errorRate=new float[sites]; 
        float[] mnErrorRate=new float[sites];
        for (int i = 0; i < errorCnt.length; i++) {
            System.out.printf("%d %d %d %d %n",i,mjCorrCnt[i],mnCorrCnt[i],errorCnt[i]);
            int tests=mjCorrCnt[i]+mnCorrCnt[i]+errorCnt[i];
            if(tests==0) {errorRate[i]=Float.NaN;}
            else {errorRate[i]=(float)errorCnt[i]/(float)tests;}
            if((mnCorrCnt[i]+errorCnt[i])==0) {mnErrorRate[i]=Float.NaN;}
            else {mnErrorRate[i]=(float)errorCnt[i]/(float)((2*mnCorrCnt[i])+errorCnt[i]);}
            System.out.printf("%d %d %d %d %d %g %g %n",i,mjCorrCnt[i],mnCorrCnt[i],errorCnt[i], tests, errorRate[i], mnErrorRate[i]);
        }
        if(!h5w.exists(HapMapHDF5Constants.IBSERRORRATE_DESC)) h5w.createFloatArray(HapMapHDF5Constants.IBSERRORRATE_DESC, errorRate.length);
        h5w.writeFloatArray(HapMapHDF5Constants.IBSERRORRATE_DESC, errorRate);
        if(!h5w.exists(HapMapHDF5Constants.IBSMINORERRORRATE_DESC)) h5w.createFloatArray(HapMapHDF5Constants.IBSMINORERRORRATE_DESC, mnErrorRate.length);
        h5w.writeFloatArray(HapMapHDF5Constants.IBSMINORERRORRATE_DESC, mnErrorRate);
    }
    
    
    private int[] addTwoVector(int[] src, int[] addendum) {
        for (int i = 0; i < addendum.length; i++) {
            src[i]+=addendum[i];
        }
        return src;
    }
    
    public static void saveAnnotationsToFile(String filename, String outFile) {
        IHDF5Reader reader = HDF5Factory.openForReading(filename);
      //  int[] variableSites = reader.readIntArray(HapMapHDF5Constants.POSITIONS);
        String delimiter="\t";
        List<HDF5LinkInformation> fields=reader.getAllGroupMemberInformation(HapMapHDF5Constants.SITE_DESC, true);
        List<HDF5LinkInformation> fields2=new ArrayList(fields);
        for (HDF5LinkInformation is : fields) {
            //if(is.isGroup()==false) continue;
            if(is.isGroup())fields2.addAll(reader.getAllGroupMemberInformation(is.getPath(), true));
        }
        float[][] fa=new float[20][];  
        String[] fNames=new String[20];
        int[][] ia=new int[20][]; 
        String[] iNames=new String[20];
        int currentFA=0;
        int currentIA=0;
        for (HDF5LinkInformation is : fields2) {
            System.out.println(is.getPath().toString()+"::"+reader.getObjectType(is.getPath()).toString());
            if(is.isDataSet()==false) continue;
            HDF5DataSetInformation info=reader.getDataSetInformation(is.getPath());
            if(info.getTypeInformation().getDataClass()==HDF5DataClass.FLOAT) {
                fNames[currentFA]=is.getName();
                fa[currentFA]=reader.readFloatArray(is.getPath());
                currentFA++;
            } else if(info.getTypeInformation().getDataClass()==HDF5DataClass.INTEGER) {
                iNames[currentIA]=is.getName();
                ia[currentIA]=reader.readIntArray(is.getPath());
                currentIA++;
            }
            
            System.out.println(is.getPath().toString()+"::"+reader.getDataSetInformation(is.getPath()).toString());
        }
        StringBuilder sb=new StringBuilder("Site"+delimiter);
        for (int fi = 0; fi < currentFA; fi++) {sb.append(fNames[fi]); sb.append(delimiter);}
        for (int ii = 0; ii < currentIA; ii++) {sb.append(iNames[ii]); sb.append(delimiter);}
        System.out.println(sb.toString());
        for (int i = 0; i < fa[0].length; i++) {
            sb=new StringBuilder();
           // sb.append(variableSites[i]);sb.append(delimiter);
            for (int fi = 0; fi < currentFA; fi++) {
                sb.append(fa[fi][i]);sb.append(delimiter);             
            }
            for (int ii = 0; ii < currentIA; ii++) {
                sb.append(ia[ii][i]);sb.append(delimiter);             
            }
            //sb.append(reader.readFloat(outFile));sb.append(delimiter);
            System.out.println(i+delimiter+sb.toString());
        }
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "FilterErrorForBiparental";
    }

    @Override
    public String getToolTipText() {
        return "Filters and estimates error rates from biparental populations";
    }
    
    private static void filterOutHets(String infile, String outfile) {
        Alignment a = ImportUtils.readFromHapmap(infile, null);
        IntArrayList keep=new IntArrayList();
        double siteD=a.getSiteCount();
        for (int i = 0; i < a.getSequenceCount(); i++) {
            double presentProp=(double)a.getTotalNotMissingForTaxon(i)/siteD;
            double hetProp=(double)a.getHeterozygousCountForTaxon(i)/siteD;
            double scaleHetProp=hetProp/presentProp;
            System.out.printf("%d %g %g %g ",i,presentProp, hetProp, scaleHetProp);
            if((presentProp>0.05)&&(scaleHetProp<0.02)) {
                keep.add(i);
                System.out.println("Keep:"+a.getFullTaxaName(i));
            } else {
                System.out.println("Delete:"+a.getFullTaxaName(i));
            }
        }
        IdGroup sid=new SimpleIdGroup(keep.size());
        for (int i = 0; i < keep.size(); i++) {
            sid.setIdentifier(i, a.getIdGroup().getIdentifier(keep.get(i)));
        }
        System.out.println("Keep size:"+keep.size());
//        System.out.println(keep.toString());
       FilterAlignment fa = (FilterAlignment)FilterAlignment.getInstance(a, sid);
       Alignment a2 = BitAlignment.getInstance(fa, true);
       ExportUtils.writeToHDF5(a2, outfile);
    }
    
    private static void filterOutHetMAFRatio(String infile, String outfile) {
        Alignment a = ImportUtils.readFromHapmap(infile, null);
        IntArrayList keep=new IntArrayList();
        double siteD=a.getSiteCount();
        int taxa=a.getSequenceCount();
        for (int i = 0; i < a.getSiteCount(); i++) {
            double maf=a.getMinorAlleleFrequency(i);
            int siteCovCnt=a.getTotalGametesNotMissing(i)/2;
            double scov=(double)siteCovCnt/(double)taxa;
            double hets=(double)a.getHeterozygousCount(i)/(double)siteCovCnt;  //terry's proportion hets is off all sites
            double hDivPresMAF=hets/(scov*maf);
            if(hDivPresMAF<0.4) {
                keep.add(i);
               // System.out.println("Keep:"+a.getFullTaxaName(i));
            } else {
               // System.out.println("Delete:"+a.getFullTaxaName(i));
            }
        }
        System.out.println("Sites Keep size:"+keep.size());
//        System.out.println(keep.toString());
       FilterAlignment fa = (FilterAlignment)FilterAlignment.getInstance(a, keep.elements());
       Alignment a2 = BitAlignment.getInstance(fa, true);
       ExportUtils.writeToHDF5(a2, outfile);
    }
    
        private static void XfilterOutHets(String infile, String outfile) {
        Alignment a = ImportUtils.readFromHapmap(infile, null);
        IntArrayList evalLowHets=new IntArrayList();
        double siteD=a.getSiteCount();
        for (int i = 0; i < a.getSequenceCount(); i++) {
            double presentProp=(double)a.getTotalNotMissingForTaxon(i)/siteD;
            double hetProp=(double)a.getHeterozygousCountForTaxon(i)/siteD;
            double scaleHetProp=hetProp/presentProp;
            System.out.printf("%d %g %g %g ",i,presentProp, hetProp, scaleHetProp);
            if((presentProp>0.05)&&(scaleHetProp<0.02)) {
                evalLowHets.add(i);
                System.out.println("Keep:"+a.getFullTaxaName(i));
            } else {
                System.out.println("Delete:"+a.getFullTaxaName(i));
            }
        }
        evalLowHets.trimToSize();
        IdGroup sid=new SimpleIdGroup(evalLowHets.size());
        for (int i = 0; i < evalLowHets.size(); i++) {
            sid.setIdentifier(i, a.getIdGroup().getIdentifier(evalLowHets.get(i)));
        }
        System.out.println("Number of Taxa Used to Count Hets:"+evalLowHets.size());
//        System.out.println(keep.toString());
       FilterAlignment fa = (FilterAlignment)FilterAlignment.getInstance(a, sid);
 //      Alignment a2 = BitAlignment.getInstance(fa, true);
 //      ExportUtils.writeToHDF5(a2, outfile);
     //   Alignment fa = ImportUtils.readFromHapmap(infile, null);
        IntArrayList keepSites=new IntArrayList();
  //      double siteD=fa.getSiteCount();
        int taxa=fa.getSequenceCount();
        for (int i = 0; i < fa.getSiteCount(); i++) {
            double maf=fa.getMinorAlleleFrequency(i);
            int siteCovCnt=fa.getTotalGametesNotMissing(i)/2;
            double scov=(double)siteCovCnt/(double)taxa;
            double hets=(double)fa.getHeterozygousCount(i)/(double)siteCovCnt;  //terry's proportion hets is off all sites
            double hDivPresMAF=hets/(scov*maf);
            if(hDivPresMAF<0.5) {
                keepSites.add(i);
               // System.out.println("Keep:"+a.getFullTaxaName(i));
            } else {
               // System.out.println("Delete:"+a.getFullTaxaName(i));
            }
        }
        keepSites.trimToSize();
        System.out.println("Sites Keep size:"+keepSites.size());
//        System.out.println(keep.toString());
       FilterAlignment fa2 = (FilterAlignment)FilterAlignment.getInstance(a, keepSites.elements());
       Alignment a2 = BitAlignment.getInstance(fa2, true);
       ExportUtils.writeToHapmap(a2, false, outfile, '\t', null);
       //ExportUtils.writeToHDF5(a2, outfile);
    }
    
    public static void testNeighborLD(String infileH5) {
        String family = "Z008";
        int window = 20;
        Alignment a = BitAlignmentHDF5.getInstance(infileH5);
        IdGroup ids = a.getIdGroup();
        double[][] avgr = new double[26][a.getSiteCount()];
       // ArrayList<Identifier> subgroup = new ArrayList<Identifier>();
        int n = ids.getIdCount();
        boolean[] keep = new boolean[n];
        for (int i = 0; i < n; i++) {
                Identifier id = ids.getIdentifier(i);
                if (id.getName().startsWith(family)) keep[i] = true;
                else keep[i] = false;
        }
        
        IdGroup subIdGroup = IdGroupUtils.idGroupSubset(ids, keep);
        Alignment b = FilterAlignment.getInstance(a, subIdGroup);
        Alignment c = BitAlignment.getInstance(b, true);
        BitSet polybits;
        polybits = NucleotideImputationUtils.whichSitesArePolymorphic(c, .8, .1);
        System.out.println("polybits cardinality = " + polybits.cardinality());

        OpenBitSet filteredBits = NucleotideImputationUtils.whichSnpsAreFromSameTag(c, polybits);
        System.out.println("filteredBits.cardinality = " + filteredBits.cardinality());

        System.out.printf("A:%d B:%d C:%d %n", a.getSiteCount(), b.getSiteCount(), c.getSiteCount());
    	int nsites = c.getSiteCount();
    	int ntests = (int) filteredBits.cardinality();
    	int[] pos = new int[ntests];
    	
    	int siteCount = 0;
    	for (int s = 0; s < nsites; s++) {
    		if (filteredBits.fastGet(s)) {
    			avgr[0][s] = NucleotideImputationUtils.neighborLD(c, s, window, filteredBits);
    		}
    	}
        for (int i = 0; i < nsites; i++) {
            System.out.printf("%d %d %g %n",i,a.getPositionInLocus(i),avgr[0][i]);
            
        }
    }
    
    private static ArrayList<ArrayList<Identifier>> classifyTaxaToPops(Alignment a, String popMask) {
        ArrayList<String> popNames = new ArrayList<String>();
        ArrayList<ArrayList<Identifier>> result=new ArrayList<ArrayList<Identifier>>();
        Pattern pat = Pattern.compile(popMask);
        IdGroup idg=a.getIdGroup();
        for (int i = 0; i < idg.getIdCount(); i++) {
            String name=idg.getIdentifier(i).getFullName();
            Matcher m = pat.matcher(name);
            if (m.find()) {
                String tn = name.substring(m.start(), m.end());
                if (!popNames.contains(tn)) {
                    popNames.add(tn);
                    result.add(new ArrayList<Identifier>());
                }
                int popOfTaxa = popNames.indexOf(tn);
                result.get(popOfTaxa).add(idg.getIdentifier(i));
            }
        }
        return result;
    }
    
    private static void calcLDByPop(String infileH5, String hdf5File) {
        Alignment a = BitAlignmentHDF5.getInstance(infileH5);
        FisherExact fe=new FisherExact(1000);
        ArrayList<ArrayList<Identifier>> pops=classifyTaxaToPops(a,"Z0[0-1][0-9]E");
        int minCntForLD = 10;
        double[][] ld34ByPop = new double[pops.size()][a.getSiteCount()];
        double[][] ldMaxByPop = new double[pops.size()][a.getSiteCount()];
        IdGroup idg = a.getIdGroup();
        for (int pop = 0; pop < pops.size(); pop++) {
            Arrays.fill(ld34ByPop[pop], Double.NaN);
            Arrays.fill(ldMaxByPop[pop], Double.NaN);
            Identifier[] ids = new Identifier[pops.get(pop).size()];
            pops.get(pop).toArray(ids);
            Alignment pa = FilterAlignment.getInstance(a, new SimpleIdGroup(ids));
            pa=BitAlignment.getInstance(pa, true);
            System.out.printf("Pop");
            IntArrayList segSites=new IntArrayList();
            for (int i = 0; i < pa.getSiteCount(); i++) {
                if((pa.getMinorAlleleCount(i)>5)&&(pa.getMajorAlleleCount(i)>20)) segSites.add(i);
            }
          //  int windowSize=10;
            int minPhysDist = 1000000;
            int numberOfTests=60;
            int sites=segSites.size();
            for (int i = 0; i < segSites.size(); i++) {
                ArrayList<Double> obsR2 = new ArrayList<Double>();
                int leftSite=i, rightSite=i, j=-1;
                int attemptTests=0, completedTests=0;
                int position=a.getPositionInLocus(segSites.get(i)), dist=-1, minSigDist=Integer.MAX_VALUE;
                BitSet rMj = pa.getAllelePresenceForAllTaxa(segSites.get(i), 0);
                BitSet rMn = pa.getAllelePresenceForAllTaxa(segSites.get(i), 1); 
//                System.out.printf("%d Site:%s Position:%d rMjCard:%d %n",i,
//                        segSites.get(i),pa.getPositionInLocus(segSites.get(i)), rMj.cardinality());
                while((completedTests<numberOfTests)&&((leftSite>0)||(rightSite+1<sites))) {
                    int rightDistance=(rightSite+1<sites)?a.getPositionInLocus(segSites.get(rightSite+1))-position:Integer.MAX_VALUE;
                    int leftDistance=(leftSite>0)?position-a.getPositionInLocus(segSites.get(leftSite-1)):Integer.MAX_VALUE;
                    if(rightDistance<leftDistance) {rightSite++; j=segSites.get(rightSite); dist=rightDistance;}
                    else {leftSite--; j=segSites.get(leftSite); dist=leftDistance;}
                    if(dist<minPhysDist) continue;
                    attemptTests++;
                    BitSet cMj = pa.getAllelePresenceForAllTaxa(j, 0);
                    BitSet cMn = pa.getAllelePresenceForAllTaxa(j, 1);
//                    System.out.printf("jw Site:%d Position:%d rMjCard:%d %n",j
//                        ,pa.getPositionInLocus(j), cMj.cardinality());
//                    LDResult result=LinkageDisequilibrium.getLDForSitePair(rMj,rMn,cMj,cMn,3,20,-1.0f,fe);
//                    if(Float.isNaN(result[1])) continue;
//                    obsR2.add((double)result[1]);
                    completedTests++;
                }
                Collections.sort(obsR2);
 //               System.out.printf("Pop:%d Site:%d LD:%s %n",pop,i,obsR2.toString());
                if (obsR2.size() > 4) {
                    ld34ByPop[pop][segSites.get(i)] = obsR2.get(obsR2.size() * 3/4);
                    ldMaxByPop[pop][segSites.get(i)] = obsR2.get(obsR2.size()-1);
//                    System.out.printf("Pop:%d Position:%d Site:%d 3/4LD:%g MaxLD:%g %n",pop,i,
//                            pa.getPositionInLocus(segSites.get(i)),ld34ByPop[pop][segSites.get(i)],obsR2.get(obsR2.size()-1));
                }
            }
            // System.out.println("POP:"+pop+Arrays.toString(ldByPop[pop]));
        }
        float[] mean34LD=new float[ldMaxByPop[0].length];
        float[] meanMaxLD=new float[ldMaxByPop[0].length];
        float[] maxMaxLD=new float[ldMaxByPop[0].length];  Arrays.fill(maxMaxLD, -1f);
        float[] minMaxLD=new float[ldMaxByPop[0].length];  Arrays.fill(minMaxLD, 2f);
        int[] cntPops=new int[ldMaxByPop[0].length];
        for (int i = 0; i < ldMaxByPop[0].length; i++) {
   //         int cntPops=0;
       //     double mean34LD=0, meanMaxLD=0, maxMaxLD=-1, minMaxLD=1.5;
            for (int j = 0; j < ldMaxByPop.length; j++) {
                if(!Double.isNaN(ld34ByPop[j][i])) {
                    cntPops[i]++;
                    mean34LD[i]+=ld34ByPop[j][i];
                    meanMaxLD[i]+=ldMaxByPop[j][i];
                    if(ldMaxByPop[j][i]>maxMaxLD[i]) maxMaxLD[i]=(float)ldMaxByPop[j][i];
                    if(ldMaxByPop[j][i]<minMaxLD[i]) minMaxLD[i]=(float)ldMaxByPop[j][i];
                }

            }
            mean34LD[i]/=(float)cntPops[i];
            meanMaxLD[i]/=(float)cntPops[i];
            if(maxMaxLD[i]<0) maxMaxLD[i]=Float.NaN;
            if(minMaxLD[i]>1) minMaxLD[i]=Float.NaN;
            System.out.printf("%d %d %d %g %g %g %g %n",i,a.getPositionInLocus(i),
                    cntPops[i],mean34LD[i],meanMaxLD[i], maxMaxLD[i], minMaxLD[i]);
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        myLogger.info("Annotating HDF5 file with BiParental: " + hdf5File);
//        config.overwrite();
//        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        if(!h5w.exists(HapMapHDF5Constants.LD_DESC)) h5w.createGroup(HapMapHDF5Constants.LD_DESC);
        if(!h5w.exists(HapMapHDF5Constants.BPLDMean34)) h5w.createFloatArray(HapMapHDF5Constants.BPLDMean34, mean34LD.length);
        h5w.writeFloatArray(HapMapHDF5Constants.BPLDMean34, mean34LD);
        if(!h5w.exists(HapMapHDF5Constants.BPLDMeanMax)) h5w.createFloatArray(HapMapHDF5Constants.BPLDMeanMax, meanMaxLD.length);
        h5w.writeFloatArray(HapMapHDF5Constants.BPLDMeanMax, meanMaxLD);
        if(!h5w.exists(HapMapHDF5Constants.BPmaxMaxLD)) h5w.createFloatArray(HapMapHDF5Constants.BPmaxMaxLD, maxMaxLD.length);
        h5w.writeFloatArray(HapMapHDF5Constants.BPmaxMaxLD, maxMaxLD);
        if(!h5w.exists(HapMapHDF5Constants.BPminMaxLD)) h5w.createFloatArray(HapMapHDF5Constants.BPminMaxLD, minMaxLD.length);
        h5w.writeFloatArray(HapMapHDF5Constants.BPminMaxLD, minMaxLD);
        if(!h5w.exists(HapMapHDF5Constants.BPPopCnt)) h5w.createIntArray(HapMapHDF5Constants.BPPopCnt, cntPops.length);
        h5w.writeIntArray(HapMapHDF5Constants.BPPopCnt, cntPops);
        h5w.close();
    }
    
    
    
//    
//     public static final String BPLDMean34 = LD_DESC+"/BPLD34MedianR2";
//    public static final String BPLDMeanMax = LD_DESC+"/BPLD34MaxR2";
//    public static final String BPmaxMaxLD = LD_DESC+"/BPmaxMaxLDR2";
//    public static final String BPminMaxLD = LD_DESC+"/BPminMaxLDR2";
    
    public static void main(String[] args) {
        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";
  //      String infile=root+"Z0NE00N_chr10S.hmp.txt.gz";
  //      String infile=root+"All_chr10S.hmp.txt.gz";
        String infile=root+"all130313.c10.hmp.txt.gz";
        String outfile=root+"all130313.c10.anno.hmp.h5";
 //       XfilterOutHets(root+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_Rev1_chr10.hmp.txt.gz",root+"noHighHSM5All_chr10.hmp.txt.gz");
//        System.exit(0);
//        filterOutHets(infile,root+"nhzAll_chr10S");
// //       filterOutHets(infile,root+"hzZ0NE00N_chr10S");
//        infile=root+"hzZ0NE00N_chr10S.hmp.txt.gz";
////        Alignment a = ImportUtils.readFromHapmap(infile, null);
////        ExportUtils.writeToHDF5(a, outfile);
//        System.exit(0);
// //       infile=outfile;
//        filterOutHetMAFRatio(infile,root+"FILTnhzAll_chr10S.hmp.txt.gz");
        saveAnnotationsToFile(outfile,root+"test.txt");
        System.exit(0);
//        Alignment align = BitAlignmentHDF5.getInstance(root+"nhzAll_chr10S.hmp.h5");
//        ExportUtils.writeToHapmap(align, false, root+"nhzAll_chr10S.hmp.txt.gz", '\t', null);
//        System.exit(0);
        
        String infileH5=root+"all130313.c10.hmp.h5";
        String annofileH5=root+"annoNoHets.hmp.h5";
 //       calcLDByPop(infileH5,outfile);
 //       System.exit(0);
 //       testNeighborLD(infileH5);
 //       System.exit(0);
        
        String[] args2 = new String[]{
            "-hmp", infile,
            "-o", outfile,
            "-mxE", "0.01",
            
//            "-sC", "10", // Start chromosome
//            "-eC", "10" // End chromosome
        };

        AnnotateGBSHDF5Plugin plugin = new AnnotateGBSHDF5Plugin();
        plugin.setParameters(args2);
        plugin.performFunction(null);
    }
}
