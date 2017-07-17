/*
 * TagsByTaxaByteFileMap
 */
package net.maizegenetics.gbs.tagdist;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;

import java.nio.ByteBuffer;

import java.util.Random;

/**
 * A byte implementation of TagByTaxa, maximum depth of tags is 127.
 * This version works directly from a RandomAccessFile, which is slow but can work
 * with very large files.  It will need to be transitioned nio2, when Java 7 comes out.
 *
 * If possible, it holds the current row cached in memory, so try to work with the files
 * sequentially.
 * @author edbuckler
 */
public class TagsByTaxaByteFileMap extends AbstractTagsByTaxa {

    private RandomAccessFile theRAF;
    private long dataStartPos = -1;  //the position of the first tag
    private long distStartPos = -1;  //the position of the first element describing the distribution of the taxa (dataStartPos+firstTagsLength+firstTagLength)
    private long bufferedTagIndex = Integer.MIN_VALUE;
    private byte[] bufferedTagDist = null;
//    private int numLongPerTaxaDist=-1;
    private long byteLenRow = -1;
    private boolean rowSetMethod = false;
    private boolean bufferChanged = false;

    public TagsByTaxaByteFileMap(String infile) {
        readDistFile(new File(infile), FilePacking.Byte);
        bufferedTagDist = new byte[taxaNum];
        getReadCountForTagTaxon(0, 0);
    }

    /**
     * Creates an empty Byte File Map
     * @param outFile
     * @param taxaNames
     * @param theTags 
     */
    public TagsByTaxaByteFileMap(String outFile, String[] taxaNames, Tags theTags) {
        int hapsOutput = 0;
        int tN = taxaNames.length;
        tagLengthInLong = theTags.getTagSizeInLong();
        byte[] obs = new byte[tN];

        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            fw.writeInt(theTags.getTagCount());
            fw.writeInt(theTags.getTagSizeInLong());
            fw.writeInt(tN);
            for (int t = 0; t < tN; t++) {
                fw.writeUTF(taxaNames[t]);
            }
            for (int i = 0; i < theTags.getTagCount(); i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(theTags.getTag(i)[j]);
                }
                fw.writeByte(theTags.getTagLength(i));
                fw.write(obs);
                hapsOutput++;
            }
            fw.close();
        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }
        readDistFile(new File(outFile), FilePacking.Byte);
        bufferedTagDist = new byte[taxaNum];
        getReadCountForTagTaxon(0, 0);
    }

    @Override
    public void readDistFile(File inFile, FilePacking binary) {
        int hapsOutput = 0;
        try {
            theRAF = new RandomAccessFile(inFile, "rw");
            int tagNum = theRAF.readInt();
            tagLengthInLong = theRAF.readInt();
            taxaNum = theRAF.readInt();
            initMatrices(taxaNum, tagNum);
            for (int t = 0; t < taxaNum; t++) {
                taxaNames[t] = theRAF.readUTF();
            }
            dataStartPos = theRAF.getFilePointer();
            System.out.printf("nBytes to end of taxa names: %d %n", dataStartPos);
            distStartPos = dataStartPos + ((long) tagLengthInLong * 8L) + 1L;  // "L" casts the number preceding it as a long (rather than an int)
            byteLenRow = 1 + taxaNum + (8 * tagLengthInLong);
            byte[] b = new byte[(int) byteLenRow];
            for (int i = 0; i < tagNum; i++) {
                if (i % 1000000 == 0) {
                    System.out.println("Read tag" + i);
                }
                theRAF.read(b);
                ByteBuffer bb = ByteBuffer.wrap(b);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = bb.getLong();
                }
                tagLength[i] = bb.get();
                hapsOutput++;
            }
        } catch (Exception e) {
            System.out.println("Catch in reading input file: " + e);
            e.printStackTrace();
        }
        System.out.println("Number of Taxa in file:" + taxaNum);
        System.out.println("Number of Haplotypes in file:" + hapsOutput);
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        if (value == getReadCountForTagTaxon(tagIndex, taxaIndex)) {
            return;
        }
        if (rowSetMethod) {
            setReadCountForTagTaxonRow(tagIndex, taxaIndex, value);
        } else {
            setReadCountForTagTaxonPoint(tagIndex, taxaIndex, value);
        }
    }

    public void setFastPointReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        setReadCountForTagTaxonPoint(tagIndex, taxaIndex, value);
    }

    public void setFastReadCountForAdjacentTagTaxonSet(int tagIndex, int firstTaxaIndex, byte[] values) {
        try {
            long pos = (byteLenRow * (long) tagIndex) + distStartPos + firstTaxaIndex;
            theRAF.seek(pos);
            theRAF.write(values);
        } catch (IOException e) {
            System.out.println("Catch in reading bufferTagDist: " + e);
            e.printStackTrace();
        }
    }

    private void setReadCountForTagTaxonRow(int tagIndex, int taxaIndex, int value) {
        if (tagIndex != bufferedTagIndex) {
            bufferTagDist(tagIndex);
        }
        if (value > Byte.MAX_VALUE) {
            bufferedTagDist[taxaIndex] = Byte.MAX_VALUE;
        } else if (value < 0) {
            bufferedTagDist[taxaIndex] = 0;
        } else {
            bufferedTagDist[taxaIndex] = (byte) value;
        }
        bufferChanged = true;
    }

    private void setReadCountForTagTaxonPoint(int tagIndex, int taxaIndex, int value) {
        try {
            long pos = (byteLenRow * (long) tagIndex) + distStartPos + taxaIndex;
            theRAF.seek(pos);
            theRAF.write(value);
        } catch (IOException e) {
            System.out.println("Catch in reading bufferTagDist: " + e);
            e.printStackTrace();
        }
    }

    public void setMethodByRows(boolean rowSetMethod) {
        this.rowSetMethod = rowSetMethod;
    }

    public void getFileReadyForClosing() {
        saveCurrentBufferTagDist();
    }

    @Override
    public int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        if (tagIndex != bufferedTagIndex) {
            bufferTagDist(tagIndex);
        }
        return bufferedTagDist[taxaIndex];
    }

    synchronized private void saveCurrentBufferTagDist() {
        try {
            //long oldpos=(byteLenRow*(long)bufferedTagIndex)+dataStartPos+((long)tagLengthInLong*8L)+1L;
            long oldpos = (byteLenRow * (long) bufferedTagIndex) + distStartPos;
            theRAF.seek(oldpos);
            theRAF.write(bufferedTagDist);
            bufferChanged = false;
        } catch (IOException e) {
            System.out.println("Catch in reading saveCurrentBufferTagDist: " + e);
            e.printStackTrace();
        }
    }

    synchronized private void bufferTagDist(int tagIndex) {
        try {
            if (bufferChanged) {//save the old buffer
                saveCurrentBufferTagDist();
            }
            //        long pos=(byteLenRow*(long)tagIndex)+dataStartPos+((long)tagLengthInLong*8L)+1L;
            long pos = (byteLenRow * (long) tagIndex) + distStartPos;
            theRAF.seek(pos);
            theRAF.read(bufferedTagDist);
            bufferedTagIndex = tagIndex;
            bufferChanged = false;
        } catch (IOException e) {
            System.out.println("Catch in reading bufferTagDist: " + e);
            e.printStackTrace();
        }
    }

    @Override
    public void initMatrices(int taxaNum, int tagNum) {
        taxaNames = new String[taxaNum];
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        throw new UnsupportedOperationException("Not supported yet.");
//        taxaNum=tagDist.length+addTaxaNames.length;
//        BitSet[] newTagDist = new BitSet[taxaNum];
//        String[] newTaxaNames = new String[taxaNum];
//        for(int i=0; i<newTagDist.length; i++) {
//            if(i<tagDist.length) {
//                newTagDist[i]=tagDist[i];
//                newTaxaNames[i]=taxaNames[i];
//                }
//            else {
//                newTagDist[i]=new BitSet(getTagCount());
//                newTaxaNames[i]=addTaxaNames[i-tagDist.length];
//            }
//        }
//        tagDist=newTagDist;
//        taxaNames=newTaxaNames;
    }

    public static void main(String[] args) {

        System.out.println("Running main method in TagsByTaxaBitFileMap");
        String infile = "/media/jvh053111/all_zea_tbt_110411/10243151_81N4HABXX_s_6.tbt.bin";
//         String infile="/Users/edbuckler/SolexaAnal/GBS/dist/NAMsort_110303.bibin";
        //    String infile="/Users/edbuckler/SolexaAnal/GBS/test/real_output.bibin";
        String outfile2 = "/media/Cache/chunk1.tbt.byte";
        TagsByTaxaByteFileMap tbtb2 = new TagsByTaxaByteFileMap(infile);
        long currTime = System.currentTimeMillis();
        int count = 0;
        TagsByTaxaByteFileMap tbtb3 = new TagsByTaxaByteFileMap(outfile2, tbtb2.getTaxaNames(), tbtb2);
        tbtb3.setMethodByRows(true);

        Random rnd = new Random();
        int rows = 10000, currRow = 10;
//        for(int i=0; i<rows; i++) {
////            if(i%4!=0) {currRow++;} else {currRow=rnd.nextInt(8000000);}
////            if(tbtb2.getReadCountForTagTaxon(currRow, i%280)==1) count++;
//
//            if(tbtb3.getReadCountForTagTaxon(75, i%280)==1) {count++;
//                System.out.println("i="+i%280);
//            }
//        }
        for (int i = 0; i < tbtb2.getTagCount(); i++) {
            for (int j = 0; j < tbtb3.getTaxaCount(); j++) {
                tbtb3.setReadCountForTagTaxon(i, j, tbtb2.getReadCountForTagTaxon(i, j));
            }
        }
        tbtb3.getFileReadyForClosing();
        tbtb2.writeDistFile(new File("/media/Cache/tbtb2.tbt.txt"), FilePacking.Text, 0);
        tbtb3.writeDistFile(new File("/media/Cache/tbtb3.tbt.txt"), FilePacking.Text, 0);
        long totTime = System.currentTimeMillis() - currTime;
        double rate = (double) rows / (double) totTime;
        System.out.printf("Total time: %d   Count:%d Rate: %g %n", totTime, count, rate);
        //      tbtb2.writeDistFile(new File(outfile2), FilePacking.Text, -1);
    }
}
