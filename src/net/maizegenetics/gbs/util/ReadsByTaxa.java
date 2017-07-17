/*
 * ReadsByTaxa
 */
package net.maizegenetics.gbs.util;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;

import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.util.BitUtil;

import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author ed
 */
public class ReadsByTaxa implements Reads {

    public long[][] haplotype;
    public byte[][] hapDist;
    public int taxaNum = 0;
    public int haplotypeNum = 0;
    public String[] taxaNames;
    public int tagLengthInLong;
    public byte[] tagLength;

    public ReadsByTaxa() {
    }

    public ReadsByTaxa(String infile, boolean binary) {
        readDistFile(new File(infile), binary);
    }

    public ReadsByTaxa(String[] taxaNames, Reads theDistinctReads) {
        this.taxaNames = taxaNames.clone();
        taxaNum = taxaNames.length;
        haplotypeNum = theDistinctReads.getReadTotal();
        haplotype = new long[2][haplotypeNum];
        for (int i = 0; i < haplotypeNum; i++) {
            long[] h = theDistinctReads.getRead(i);
            haplotype[0][i] = h[0];
            haplotype[1][i] = h[1];
        }
        hapDist = new byte[haplotypeNum][taxaNum];
    }

    public ReadsByTaxa(long[][] reads, byte[][] readDist, String[] namesForTaxa) {
        haplotype = reads;
        hapDist = readDist;
        taxaNames = namesForTaxa;
        taxaNum = namesForTaxa.length;
        haplotypeNum = reads[0].length;
    }

    public int getTaxaCount() {
        return taxaNames.length;
    }

    public String getTaxaName(int taxaIndex) {
        return taxaNames[taxaIndex];
    }

    public String[] getTaxaNames() {
        return taxaNames;
    }

    public int getIndexOfTaxaName(String taxon) {
        for (int i = 0; i < taxaNames.length; i++) {
            if (taxon.equals(taxaNames[i])) {
                return i;
            }
        }
        return -1;
    }

    public int getReadCountForTaxa(int readIndex, int taxaIndex) {
        return (int) hapDist[readIndex][taxaIndex];
    }

    public void setReadCountForTaxa(int readIndex, int taxaIndex, int value) {
        if (value > Byte.MAX_VALUE) {
            hapDist[readIndex][taxaIndex] = Byte.MAX_VALUE;
        } else if (value < 0) {
            hapDist[readIndex][taxaIndex] = 0;
        } else {
            hapDist[readIndex][taxaIndex] = (byte) value;
        }
    }

    public void addToReadCountForTaxa(int readIndex, int taxaIndex, int addValue) {
        setReadCountForTaxa(readIndex, taxaIndex, addValue + hapDist[readIndex][taxaIndex]);
    }

    public byte[] getReadCountsForTaxa(int readIndex) {
        return hapDist[readIndex].clone();
    }

    public int getTaxaCountForRead(int readIndex) {   // how many taxa was a given read seen in?
        int nTaxaWData = 0;
        for (int cnt : hapDist[readIndex]) {
            if (cnt > 0) {
                ++nTaxaWData;
            }
        }
        return nTaxaWData;
    }

    @Override
    public long[] getRead(int i) {
        if (i >= haplotypeNum) {
            return null;
        }
        long[] result = {haplotype[0][i], haplotype[1][i]};
        return result;
    }

    /**
     * Gets the first index of a read (the only one if a unique list).
     * If the read is not found then it return
     * a negative value indicating its insertion point.
     * @param read as a compressed long array
     * @return the index of the read in the array
     */
    @Override
    public int getReadIndex(long[] read) {
        int hit = Arrays.binarySearch(haplotype[0], read[0]);
        if (hit < 1) {
            return hit;
        }
        while (haplotype[0][hit - 1] == read[0]) {
            hit--;
        }
        while ((haplotype[0][hit] == read[0]) && (hit < haplotype[0].length - 1) && (haplotype[1][hit] < read[1])) {
            hit++;
        }
        if (((haplotype[0][hit] == read[0]) && (haplotype[1][hit] == read[1]))) {
            return hit;
        }
        return -hit;
    }

    /**Converts a TBT Bit file to a ReadsByTaxa object (for compatibility only).*/
    public void readTBTFile(File inFile) {
        System.out.println("Reading Haplotypes distribution from:" + inFile.toString());
        int hapsOutput = 0;
        try {
            DataInputStream rw = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 4000000));
            haplotypeNum = rw.readInt();
            tagLengthInLong = rw.readInt();
            taxaNum = rw.readInt();
            taxaNames = new String[taxaNum];

            //Init matrices
            haplotype = new long[2][haplotypeNum];
            hapDist = new byte[haplotypeNum][taxaNum];
            tagLength = new byte[haplotypeNum];

            for (int t = 0; t < taxaNum; t++) {
                taxaNames[t] = rw.readUTF();
            }

            int numberOfLongs = BitUtil.bits2words(taxaNum);
            long[] distInLong = new long[numberOfLongs];
            OpenBitSet obs;
            obs = new OpenBitSet(distInLong, taxaNum);

            for (int i = 0; i < haplotypeNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    haplotype[j][i] = rw.readLong();
                }
                tagLength[i] = rw.readByte();
                for (int j = 0; j < numberOfLongs; j++) {
                    distInLong[j] = rw.readLong();
                }

                for (int t = 0; t < taxaNum; t++) {
                    if (obs.fastGet(t)) {
                        hapDist[i][t] = 1;
                    } else {
                        hapDist[i][t] = 0;
                    }
                }
                hapsOutput++;
            }

            rw.close();
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    void readDistFile(File inFile, boolean binary) {
        System.out.println("Reading Haplotypes distribution from:" + inFile.toString());
        int hapsOutput = 0;
        if (binary) {
            try {
                DataInputStream rw = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 4000000));
                taxaNum = rw.readInt();
                haplotypeNum = rw.readInt();
                taxaNames = new String[taxaNum];
                haplotype = new long[2][haplotypeNum];
                hapDist = new byte[haplotypeNum][taxaNum];
                for (int t = 0; t < taxaNum; t++) {
                    taxaNames[t] = rw.readUTF();
                }
                for (int i = 0; i < haplotypeNum; i++) {
                    haplotype[0][i] = rw.readLong();
                    haplotype[1][i] = rw.readLong();
                    for (int t = 0; t < taxaNum; t++) {
                        hapDist[i][t] = rw.readByte();
                    }
                    hapsOutput++;
                }
                rw.close();
            } catch (Exception e) {
                System.out.println("Catch in writing output file e=" + e);
            }
        } else {
            ArrayList<String> inputLine;
            long[] haplo = new long[2];
            try {
                BufferedReader br = new BufferedReader(new FileReader(inFile), 65536);
                inputLine = new ArrayList<String>(Arrays.asList(br.readLine().split("\t")));
                taxaNum = Integer.parseInt(inputLine.get(0));
                haplotypeNum = Integer.parseInt(inputLine.get(1));
                taxaNames = new String[taxaNum];
                haplotype = new long[2][haplotypeNum];
                hapDist = new byte[haplotypeNum][taxaNum];
                inputLine = new ArrayList<String>(Arrays.asList(br.readLine().split("\t")));
                for (int t = 0; t < taxaNum; t++) {
                    taxaNames[t] = inputLine.get(t + 1);  // blank cell before the taxa list
                }
                for (int i = 0; i < haplotypeNum; i++) {
                    inputLine = new ArrayList<String>(Arrays.asList(br.readLine().split("\t")));
                    haplo = BaseEncoder.getLongArrayFromSeq(inputLine.get(0));
                    haplotype[0][i] = haplo[0];
                    haplotype[1][i] = haplo[1];
                    for (int t = 0; t < taxaNum; t++) {
                        hapDist[i][t] = Byte.valueOf(inputLine.get(t + 1));
                    }
                    hapsOutput++;
                }
            } catch (Exception e) {
                System.out.println("Catch in writing output file e=" + e);
            }
        }
        System.out.println("Number of Taxa in file:" + taxaNum);
        System.out.println("Number of Haplotypes in file:" + hapsOutput);
    }

    public void writeDistFile(File outFile, boolean binary, int minCount) {
        int hapsOutput = 0;
        int outReads = readsWCountsGreaterThanMin(minCount);
        System.out.println(outReads + " reads will be output to " + outFile.getName());
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            if (binary) {
                fw.writeInt(taxaNum);
                fw.writeInt(outReads);
                for (int t = 0; t < taxaNum; t++) {
                    fw.writeUTF(taxaNames[t]);
                }
            } else {
                fw.writeBytes(taxaNum + "\t" + outReads + "\n");
                for (int t = 0; t < taxaNum; t++) {
                    fw.writeBytes("\t" + taxaNames[t]);
                }
                fw.writeBytes("\n");
            }
            for (int i = 0; i < haplotype[0].length; i++) {
                if (getReadCount(i) < minCount) {
                    continue;
                }
                if (!binary) {
                    fw.writeBytes(
                            BaseEncoder.getSequenceFromLong(haplotype[0][i])
                            + BaseEncoder.getSequenceFromLong(haplotype[1][i]) + "\t");
                    for (int t = 0; t < taxaNum; t++) {
                        fw.writeBytes(hapDist[i][t] + "\t");
                    }
                    fw.writeBytes("\n");
                } else {
                    fw.writeLong(haplotype[0][i]);
                    fw.writeLong(haplotype[1][i]);
                    for (int t = 0; t < taxaNum; t++) {
                        fw.writeByte(hapDist[i][t]);
                    }
                }
                hapsOutput++;
            }
            fw.flush();
            fw.close();
            System.out.println("Haplotypes written to:" + outFile.toString());
            System.out.println("Number of Haplotypes in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    public void filterForListOfReads(File readsToKeepFile, File outFile, boolean binary) {
        int hapsOutput = 0;
        String readToKeep;
        ArrayList<String> readsToKeepArrayList = new ArrayList<String>();
        int outReads = 0;

        try {
            BufferedReader br = new BufferedReader(new FileReader(readsToKeepFile), 65536);
            while ((readToKeep = br.readLine()) != null) {
                readsToKeepArrayList.add(readToKeep);
            }
            String[] readsToKeep = readsToKeepArrayList.toArray(new String[readsToKeepArrayList.size()]);
            Arrays.sort(readsToKeep);

            for (int i = 0; i < haplotype[0].length; i++) {
                String readStr = BaseEncoder.getSequenceFromLong(haplotype[0][i]) + BaseEncoder.getSequenceFromLong(haplotype[1][i]);
                if (Arrays.binarySearch(readsToKeep, readStr) > -1) {
                    ++outReads;
                }
            }

            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            if (binary) {
                fw.writeInt(taxaNum);
                fw.writeInt(outReads);
                for (int t = 0; t < taxaNum; t++) {
                    fw.writeUTF(taxaNames[t]);
                }
            } else {
                fw.writeBytes(taxaNum + "\t" + outReads + "\n");
                for (int t = 0; t < taxaNum; t++) {
                    fw.writeBytes("\t" + taxaNames[t]);
                }
                fw.writeBytes("\n");
            }
            for (int i = 0; i < haplotype[0].length; i++) {
                String readStr = BaseEncoder.getSequenceFromLong(haplotype[0][i]) + BaseEncoder.getSequenceFromLong(haplotype[1][i]);
                if (Arrays.binarySearch(readsToKeep, readStr) < 0) {
                    continue;
                }
                if (!binary) {
                    fw.writeBytes(
                            BaseEncoder.getSequenceFromLong(haplotype[0][i])
                            + BaseEncoder.getSequenceFromLong(haplotype[1][i]) + "\t");
                    for (int t = 0; t < taxaNum; t++) {
                        fw.writeBytes(hapDist[i][t] + "\t");
                    }
                    fw.writeBytes("\n");
                } else {
                    fw.writeLong(haplotype[0][i]);
                    fw.writeLong(haplotype[1][i]);
                    for (int t = 0; t < taxaNum; t++) {
                        fw.writeByte(hapDist[i][t]);
                    }
                }
                hapsOutput++;
            }
            fw.flush();
            fw.close();
            System.out.println("Haplotypes written to:" + outFile.toString());
            System.out.println("Number of Haplotypes in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing filtered output file e=" + e);
        }
    }

    protected void writeReadCountFile(File outFile, boolean binary, int minCount) {
        int hapsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 4000000));
            for (int i = 0; i < haplotype[0].length - 1; i++) {
                if (getReadCount(i) >= minCount) {
                    if (!binary) {
                        fw.writeBytes(
                                BaseEncoder.getSequenceFromLong(haplotype[0][i])
                                + BaseEncoder.getSequenceFromLong(haplotype[1][i]) + " " + getReadCount(i) + "\n");
                    } else {
                        fw.writeLong(haplotype[0][i]);
                        fw.writeLong(haplotype[1][i]);
                        fw.writeInt(getReadCount(i));
                    }
                    hapsOutput++;
                }
            }
            fw.flush();
            fw.close();
            System.out.println("Reads written to:" + outFile.toString());
            System.out.println("Number of Reads in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    private int readsWCountsGreaterThanMin(int minCount) {
        int sum = 0;
        for (int i = 0; i < getReadTotal(); i++) {
            if (getReadCount(i) >= minCount) {
                sum++;
            }
        }
        return sum;
    }

    @Override
    public int getReadCount(int index) {
        int sum = 0;
        for (int cnt : hapDist[index]) {
            sum += cnt;
        }
        return sum;
    }

    @Override
    public int[] getReadIndexSet(long[] read) {
        int r = getReadIndex(read);
        if (r < 0) {
            return null;
        }
        int[] result = {r};
        return result;
    }

    @Override
    public boolean areReadsUnique() {
        return true;
    }

    @Override
    public int getReadTotal() {
        return haplotype[0].length;
    }
}
