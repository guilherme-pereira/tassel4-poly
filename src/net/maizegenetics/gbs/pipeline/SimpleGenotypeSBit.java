/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 * Store sBit for tag genetic mapping and LD detection
 * @author Fei Lu
 */
public class SimpleGenotypeSBit {
    public int taxaNum;
    public int chrNum;
    public int siteNum;
    public int wordNum;
    
    public String[] taxaNames;
    public int[] chromosomeNumber;
    /**First SNP index of a chromosome in the whole SNP list */
    public int[] chrStartIndex;
    /**Last SNP index of a chromosome in the whole SNP list, exclusive */
    public int[] chrEndIndex;
    
    public int[] position;
    public double[] maf;
    public OpenBitSet[] obsMajor;
    public OpenBitSet[] obsMinor;
    
    /**HDF5 setting*/
    private static final int BITS_TO_SHIFT_FOR_CHUNK = 8;
    private static final int CHUNK_SIZE = 1 << BITS_TO_SHIFT_FOR_CHUNK;
    int maxTaxaNameLength = 100;
    /**Attributes*/
    String ROOT = "/";
    String TAXANUM = "taxaNum";
    String CHRNUM = "chrNum";
    String SITENUM = "siteNum";
    String WORDNUM = "wordNum";
    /**Path*/
    String TAXANAMES = "taxaNames";
    String CHROMOSOME = "chromosome";
    String CHRSTARTINDEX = "chrStartIndex";
    String CHRENDINDEX = "chrEndIndex";
    String POSITION = "position";
    String MAF = "maf";
    String OBSMAJOR = "obsMajor";
    String OBSMINOR = "obsMinor";
    /**compression level*/
    HDF5IntStorageFeatures intFeature = HDF5IntStorageFeatures.createDeflation(9);
    HDF5GenericStorageFeatures genericFeature = HDF5GenericStorageFeatures.createDeflation(9);
    HDF5FloatStorageFeatures floatFeature = HDF5FloatStorageFeatures.createDeflation(9);
    
    /**
     * Convert HDF5 Alignment/Genotype file to SimpleGenotypeSBit
     * @param genotypeH5FileS
     * @param sBitFileS 
     */
    public SimpleGenotypeSBit (String genotypeH5FileS, String sBitFileS) {
        long lastTimePoint  = System.nanoTime();
        Alignment a = ImportUtils.readGuessFormat(genotypeH5FileS, true);
        taxaNames = new String[a.getTaxaCount()];
        for (int i = 0; i < taxaNames.length; i++) {
            taxaNames[i] = a.getFullTaxaName(i);
        }
        taxaNum = taxaNames.length;
        int[] chrOffSet = a.getLociOffsets();
        chromosomeNumber = new int[a.getLoci().length];
        chrStartIndex = new int[chromosomeNumber.length];
        chrEndIndex = new int[chromosomeNumber.length];
        for (int i = 0; i < chromosomeNumber.length; i++) {
            chromosomeNumber[i] = a.getLoci()[i].getChromosomeNumber();
            chrStartIndex[i] = chrOffSet[i];
            chrEndIndex[i] = chrOffSet[i] + a.getLocusSiteCount(a.getLoci()[i]);
        }
        chrNum = this.chromosomeNumber.length;
        position = a.getPhysicalPositions();
        siteNum = this.position.length;
        System.out.println("This genotype has " + taxaNum + " taxa, " + chromosomeNumber.length + " chromosomes, " + siteNum + " sites");
        System.out.println("Will be transformed to sBit with " + this.getChunkNum() + " chunks, each chunk has " + this.getChunkSize()+" sites");
        maf = new double[a.getSiteCount()];
        obsMajor = new OpenBitSet[a.getSiteCount()];
        obsMinor = new OpenBitSet[a.getSiteCount()];    
        OpenBitSet bsMa, bsMi, het;
        long majorCount, minorCount;
        a = AlignmentUtils.optimizeForSites(a);
        for (int i = 0; i < position.length; i++) {
            het = new OpenBitSet (a.getAllelePresenceForAllTaxa(i, 0).getBits().clone());
            bsMa = new OpenBitSet (a.getAllelePresenceForAllTaxa(i, 0).getBits());
            bsMi = new OpenBitSet (a.getAllelePresenceForAllTaxa(i, 1).getBits());
            het.and(bsMi);
            bsMa.xor(het);
            bsMi.xor(het);
            obsMajor[i] = bsMa;
            obsMinor[i] = bsMi;
            majorCount = bsMa.cardinality();
            minorCount = bsMi.cardinality();
            maf[i] = (double)minorCount/(majorCount+minorCount);
        }
        wordNum = obsMajor[0].getNumWords();
        System.out.println("Transform Genotype to sBit took " + this.getTimeSpanSecond(lastTimePoint) +  " seconds");
        this.writeH5File(sBitFileS);
    }
    
    /**
     * Read in HDF5 SimpleGenotypeSBit
     * @param inputfileS 
     */
    public SimpleGenotypeSBit (String inputfileS) {
        this.readH5File(inputfileS);
    }
    
    public void writeH5File (String outputFileS) {
        long lastTimePoint  = this.getCurrentTimeNano();
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(outputFileS));
        config.overwrite();
        config.useUTF8CharacterEncoding();
        IHDF5Writer h5 = config.writer();
        h5.setIntAttribute(ROOT, TAXANUM, taxaNum);
        h5.setIntAttribute(ROOT, CHRNUM, chrNum);
        h5.setIntAttribute(ROOT, SITENUM, siteNum);
        h5.setIntAttribute(ROOT, WORDNUM, wordNum);
        h5.createStringArray(TAXANAMES, maxTaxaNameLength, taxaNames.length, genericFeature);
        h5.writeStringArray(TAXANAMES, taxaNames, genericFeature);
        h5.createIntArray(CHROMOSOME, chrNum, intFeature);
        h5.writeIntArray(CHROMOSOME, chromosomeNumber, intFeature);
        h5.createIntArray(CHRSTARTINDEX, chrNum, intFeature);
        h5.writeIntArray(CHRSTARTINDEX, chrStartIndex, intFeature);
        h5.createIntArray(CHRENDINDEX, chrNum, intFeature);
        h5.writeIntArray(CHRENDINDEX, chrEndIndex, intFeature);
        h5.createIntArray(POSITION, siteNum, intFeature);
        h5.writeIntArray(POSITION, position, intFeature);
        h5.createDoubleArray(MAF, siteNum, floatFeature);
        h5.writeDoubleArray(MAF, maf, floatFeature);
        long[][] dis = new long[siteNum][];
        for (int i = 0; i < siteNum; i++) {
            dis[i] = obsMajor[i].getBits();
        }
        h5.createLongMatrix(OBSMAJOR, siteNum, wordNum, this.getChunkSize(), wordNum, intFeature);
        for (int i = 0; i < this.getChunkNum(); i++) {
            long[][] chunk = this.getSubLongMatrix(dis, this.getChunkSize(), i);
            h5.writeLongMatrixBlock(OBSMAJOR, chunk, i, 0);
        }
        for (int i = 0; i < siteNum; i++) {
            dis[i] = obsMinor[i].getBits();
        }
        h5.createLongMatrix(OBSMINOR, siteNum, wordNum, this.getChunkSize(), wordNum, intFeature);
        for (int i = 0; i < this.getChunkNum(); i++) {
            long[][] chunk = this.getSubLongMatrix(dis, this.getChunkSize(), i);
            h5.writeLongMatrixBlock(OBSMINOR, chunk, i, 0);
        }
        h5.flush();
        h5.close();
        System.out.println("Write to SimpleGenotypeSBit took " + this.getTimeSpanSecond(lastTimePoint) +  " seconds");
    }
    
    private long[][] getSubLongMatrix (long[][] dis, int actualChunkSize, int blockNumberX) {
        long[][] result = new long[actualChunkSize][];
        int chunkStartSiteIndex = actualChunkSize*blockNumberX;
        if (chunkStartSiteIndex+actualChunkSize > this.siteNum) {
            actualChunkSize = -(chunkStartSiteIndex-this.siteNum);
            for (int i = actualChunkSize; i < this.getChunkSize(); i++) {
                result[i] = new long[dis[0].length];
                for (int j = 0; j < dis[0].length; j++) {
                    result[i][j] = Long.MIN_VALUE;
                }
            }
        }
        for (int i = 0; i < actualChunkSize; i++) {
            result[i] = dis[chunkStartSiteIndex+i];
        }
        return result;
    }
    
    
    private long[][] readInSBitChunk (IHDF5Writer h5, String path) {
        long[][] dis = new long[this.siteNum][this.wordNum];
        long[][] sub = new long[this.getChunkSize()][this.wordNum];
        int actualSize = this.getChunkSize();
        int chunkStartSiteIndex;
        for (int i = 0; i < this.getChunkNum(); i++) {
            chunkStartSiteIndex = this.getChunkSize()*i;
            actualSize = this.getChunkSize();
            sub = h5.readLongMatrixBlock(path, this.getChunkSize(), this.wordNum, i, 0);
            if (chunkStartSiteIndex + this.getChunkSize() > this.siteNum) {
                actualSize = siteNum-chunkStartSiteIndex;
            }
            for (int j = 0; j < actualSize; j++) {
                dis[chunkStartSiteIndex+j] = sub[j];
            }
        }
        return dis;
    }
    
    private class MTReadInSBitChunk implements Runnable {
        IHDF5Writer h5;
        String path;
        long[][] dis;
        int chunkIndex;
        public MTReadInSBitChunk (IHDF5Writer h5, String path, long[][] dis, int chunkIndex) {
            this.h5 = h5;
            this.path = path;
            this.dis = dis;
            this.chunkIndex = chunkIndex;
        }
        
        @Override
        public void run() {
            int chunkStartSiteIndex = getChunkSize()*chunkIndex;
            int actualSize = getChunkSize();
            long[][] sub = h5.readLongMatrixBlock(path, getChunkSize(), wordNum, chunkIndex, 0);
            if (chunkStartSiteIndex + getChunkSize() > siteNum) {
                actualSize = siteNum-chunkStartSiteIndex;
            }
            for (int i = 0; i < actualSize; i++) {
                dis[chunkStartSiteIndex+i] = sub[i];
            }
            h5.close();
        }
        
    }
    
    public void readH5File (String inputFileS) {
        long lastTimePoint = this.getCurrentTimeNano();
        IHDF5Writer h5 = HDF5Factory.open(inputFileS);
        taxaNum = h5.getIntAttribute(ROOT, TAXANUM);
        chrNum = h5.getIntAttribute(ROOT, CHRNUM);
        siteNum = h5.getIntAttribute(ROOT, SITENUM);
        wordNum = h5.getIntAttribute(ROOT, WORDNUM);
        taxaNames = h5.readStringArray(TAXANAMES);
        chromosomeNumber = h5.readIntArray(CHROMOSOME);
        chrStartIndex = h5.readIntArray(CHRSTARTINDEX);
        chrEndIndex = h5.readIntArray(CHRENDINDEX);
        position = h5.readIntArray(POSITION);
        maf = h5.readDoubleArray(MAF);
        obsMajor = new OpenBitSet[siteNum];
        obsMinor = new OpenBitSet[siteNum];

        
//MT readin##############################################/        
        h5.close();
        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        long[][] dis = new long[siteNum][];
        for (int i = 0; i < this.getChunkNum(); i++) {
            pool.execute(new MTReadInSBitChunk(HDF5Factory.open(inputFileS), OBSMAJOR, dis, i));
            //new MTReadInSBitChunk(h5, OBSMAJOR, dis, i).run();
        }
        pool.shutdown();
        try {
            pool.awaitTermination(Integer.MAX_VALUE, TimeUnit.SECONDS);
        } 
        catch (InterruptedException e) {
           System.out.println(e.toString());
        }
        for (int i = 0; i < siteNum; i++) {
            obsMajor[i] = new OpenBitSet(dis[i]);
        }
        
        pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        dis = new long[siteNum][];
        for (int i = 0; i < this.getChunkNum(); i++) {
            pool.execute(new MTReadInSBitChunk(HDF5Factory.open(inputFileS), OBSMINOR, dis, i));
            //new MTReadInSBitChunk(h5, OBSMINOR, dis, i).run();
        }
        pool.shutdown();
        try {
            pool.awaitTermination(Integer.MAX_VALUE, TimeUnit.SECONDS);
        } 
        catch (InterruptedException e) {
           System.out.println(e.toString());
        }
        for (int i = 0; i < siteNum; i++) {
            obsMinor[i] = new OpenBitSet(dis[i]);
        }
//#################################################################/
/*        
//Readin###########################################################/        
        long[][] dis = this.readInSBitChunk(h5, OBSMAJOR);
        for (int i = 0; i < siteNum; i++) {
            obsMajor[i] = new OpenBitSet(dis[i]);
        }
        dis = this.readInSBitChunk(h5, OBSMINOR);
        for (int i = 0; i < siteNum; i++) {
            obsMinor[i] = new OpenBitSet(dis[i]);
        }
        h5.close();
//#################################################################/
*/ 
        System.out.println("Read SimpleGenotypeSBit took " + this.getTimeSpanSecond(lastTimePoint) +  " seconds");
    }
    
    public int getSiteNum () {
        return this.siteNum;
    }
    
    public int getTaxaNum () {
        return this.taxaNum;
    }
    
    public int getTaxonIndex (String taxonName) {
        for (int i = 0; i < taxaNum; i++) {
            if (taxaNames[i].equals(taxonName)) return i;
        }
        return -1;
    }
    
    private double getTimeSpanSecond (long lastTimePoint) {
        return (double)this.getTimeSpanNano(lastTimePoint)/1000000000;
    }
    
    private long getTimeSpanNano (long lastTimePoint) {
        return this.getCurrentTimeNano()- lastTimePoint;
    }
    
    private long getCurrentTimeNano () {
        return System.nanoTime();
    }
    
    /**
     * Return index of nearest site of given position
     * @param chr
     * @param pos
     * @return 
     */
    public int getSiteIndex (int chr, int pos) {
        int chrIndex = Arrays.binarySearch(this.chromosomeNumber, chr);
        if (chrIndex < 0) return Integer.MIN_VALUE;
        int siteIndex = Arrays.binarySearch(position, chrStartIndex[chrIndex], chrEndIndex[chrIndex], pos);
        if (siteIndex < 0) {
            siteIndex = - siteIndex - 2;
        }
        if (siteIndex < chrStartIndex[chrIndex]) siteIndex = chrStartIndex[chrIndex];
        if (siteIndex > chrEndIndex[chrIndex]) siteIndex = chrEndIndex[chrIndex];
        return siteIndex;
    }
    
    public int[] getAdjacentSiteIndexRange (int siteIndex, int siteNum) {
        int currentChr = this.getChr(siteIndex);
        int half = siteNum/2;
        int[] indexRange = new int[2];
        indexRange[0] = siteIndex-half;
        indexRange[1] = siteIndex+half+1;
        if (indexRange[0] < 0 || this.getChr(indexRange[0]) != currentChr) {
            indexRange[0] = this.chrStartIndex[this.getChrIndex(currentChr)];
        }
        if (indexRange[1] >= this.getSiteNum() || this.getChr(indexRange[1]) != currentChr) {
            indexRange[1] = this.chrEndIndex[this.getChrIndex(currentChr)];
        }
        return indexRange;
    }
    
    public int getChrIndex (int chr) {
        return Arrays.binarySearch(this.chromosomeNumber, chr);
    }
    
    public int getChr (int index) {
        int chrIndex = Arrays.binarySearch(this.chrStartIndex, index);
        if (chrIndex < 0) chrIndex = -chrIndex - 2;
        return this.chromosomeNumber[chrIndex];
    }
    
   /**
     * Return the total number of chunks
     * @return 
     */
    public int getChunkNum () {
        int num = siteNum/this.getChunkSize();
        if (siteNum% this.getChunkSize() == 0) return num;
        else return num+1;
    }
    
    /**
     * Return the chunk size (Number of tags in a chunk)
     * @return 
     */
    public int getChunkSize () {
        return this.CHUNK_SIZE;
    }
    
    public int getPosition (int siteIndex) {
        return position[siteIndex];
    }
    
    public int getChromosomeNumber (int siteIndex) {
        int hit = Arrays.binarySearch(this.chrStartIndex, siteIndex);
        if (hit < 0) hit = -hit-2;
        return this.chromosomeNumber[hit];
    }
    
    public void writeBinaryFile (String outputFileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileS), 65536));
            dos.writeInt(taxaNames.length);
            dos.writeInt(chrNum);
            dos.writeInt(siteNum);
            dos.write(wordNum);
            for (int i = 0; i < taxaNames.length; i++) {
                dos.writeUTF(taxaNames[i]);
            }
            for (int i = 0; i < this.chromosomeNumber.length; i++) {
                dos.writeInt(this.chromosomeNumber[i]);
                dos.writeInt(this.chrStartIndex[i]);
                dos.writeInt(this.chrEndIndex[i]);
            }
            long[] bits;
            for (int i = 0; i < maf.length; i++) {
                dos.writeInt(position[i]);
                dos.writeDouble(maf[i]);
                bits = obsMajor[i].getBits();
                for (int j = 0; j < bits.length; j++) {
                    dos.writeLong(bits[j]);
                }
                bits = obsMinor[i].getBits();
                for (int j = 0; j < bits.length; j++) {
                    dos.writeLong(bits[j]);
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void readBinaryFile (String inputFileS) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputFileS), 65536));
            taxaNum = dis.readInt();
            chrNum = dis.readInt();
            siteNum = dis.readInt();
            wordNum = dis.readInt();
            this.initialize();
            for (int i = 0; i < taxaNum; i++) {
                taxaNames[i] = dis.readUTF();
            }
            for (int i = 0; i < chrNum; i++) {
                this.chromosomeNumber[i] = dis.readInt();
                this.chrStartIndex[i] = dis.readInt();
                this.chrEndIndex[i] = dis.readInt();
            }
            long[] bits = new long[wordNum];
            for (int i = 0; i < siteNum; i++) {
                position[i] = dis.readInt();
                maf[i] = dis.readDouble();
                for (int j = 0; j < wordNum; j++) {
                    bits[j] = dis.readLong();
                }
                obsMajor[i] = new OpenBitSet(bits);
                for (int j = 0; j < wordNum; j++) {
                    bits[j] = dis.readLong();
                }
                obsMinor[i] = new OpenBitSet(bits);
            }
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    private void initialize () {
        taxaNames = new String[taxaNum];
        chromosomeNumber = new int[chrNum];
        chrStartIndex = new int[chrNum];
        chrEndIndex = new int[chrNum];
        position = new int[siteNum];
        maf = new double[siteNum];
        obsMajor = new OpenBitSet[siteNum];
        obsMinor = new OpenBitSet[siteNum];
    }
}
