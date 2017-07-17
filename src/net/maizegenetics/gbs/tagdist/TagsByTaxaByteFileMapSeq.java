/*
 * TagsByTaxaByteFileMapSeq
 */
package net.maizegenetics.gbs.tagdist;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import java.util.Arrays;

/**
 * A bit implementation of TagByTaxa, so only presence or absence of a tags is known.
 * This version works directly from a RandomAccessFile, which is slow but can work
 * with very large files.  It will need to be transitioned nio2, when Java 7 comes out.
 *
 * If possible, it holds the current row cached in memory, so try to work with the files
 * sequentially.
 * @author edbuckler
 */
public class TagsByTaxaByteFileMapSeq extends AbstractTagsByTaxa {

    private String name;
    private DataInputStream theRAF;
    private long bufferedTagIndex = Integer.MIN_VALUE;
    private byte[] bufferedTagDist = null;
    private long[] currTag;
    private byte currTagLength;
//    private int numLongPerTaxaDist=-1;
    private long byteLenRow = -1;
    String infileName = "Unknown";

    public TagsByTaxaByteFileMapSeq(String infile) {
        this.infileName = infile;
        readDistFile(new File(infile), FilePacking.Byte);
        bufferedTagIndex = -1;
        bufferedTagDist = new byte[taxaNum];
        currTag = new long[tagLengthInLong];
        Arrays.fill(currTag, Long.MIN_VALUE);
        name = infile;
        //  getReadCountForTagTaxon(0, 0);
    }

    @Override
    public void readDistFile(File inFile, FilePacking binary) {
        int hapsOutput = 0;
        try {
            theRAF = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 4000000));
            int tagNum = theRAF.readInt();
            tagLengthInLong = theRAF.readInt();
            taxaNum = theRAF.readInt();
            initMatrices(taxaNum, tagNum);
            for (int t = 0; t < taxaNum; t++) {
                taxaNames[t] = theRAF.readUTF();
            }
            hapsOutput = tagNum;
        } catch (Exception e) {
            System.out.println("Catch in reading input file: " + e);
            e.printStackTrace();
        }
        System.out.println("Number of Taxa in file:" + taxaNum);
        System.out.println("Number of Haplotypes in file:" + hapsOutput);
    }

    @Override
    public int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        if (tagIndex == bufferedTagIndex) {
            return bufferedTagDist[taxaIndex];
        }
        throw new UnsupportedOperationException("Not supported yet.");
    }

    synchronized public byte[] advanceToTagDist(long[] tag) {
        if (theRAF == null) {
            return (new byte[taxaNum]);
        }
        try {
            while (compareTags(currTag, tag) < 0) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    currTag[j] = theRAF.readLong();
                }
                currTagLength = theRAF.readByte();
                bufferedTagDist = new byte[taxaNum];
                theRAF.read(bufferedTagDist);
                bufferedTagIndex++;
            }
            if (compareTags(currTag, tag) == 0) {
                return bufferedTagDist;
            }
            return (new byte[taxaNum]);
        } catch (IOException e) {
            System.out.println("Closing: " + infileName);
            try {
                theRAF.close();
                theRAF = null;
            } catch (IOException ee) {
                System.out.println("Uncloseable:" + infileName);
                e.printStackTrace();
            }
            return (new byte[taxaNum]);
        }
    }

    @Override
    public void initMatrices(int taxaNum, int tagNum) {
        taxaNames = new String[taxaNum];
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setMethodByRows(boolean rowSetMethod) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void getFileReadyForClosing() {
        if (theRAF == null) {
            return;
        }
        try {
            theRAF.close();
        } catch (IOException e) {
            System.out.println("Error closing: " + infileName + e);
            e.printStackTrace();
        }
    }
}
