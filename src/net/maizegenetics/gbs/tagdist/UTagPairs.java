/*
 * UTagPairs
 */
package net.maizegenetics.gbs.tagdist;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 *
 * @author Fei Lu
 */
public class UTagPairs {

    int tagLengthInLong;
    long[][] tags;
    byte[] tagLength;
    int[] order;

    public UTagPairs(String infileS) {
        this.readTagPair(infileS);
    }

    public int getTagNum() {
        return tags[0].length;
    }

    public byte getTagLength(int index) {
        return tagLength[index];
    }

    public long[] getTag(int index) {
        long[] tag = new long[tagLengthInLong];
        for (int i = 0; i < tag.length; i++) {
            tag[i] = tags[i][index];
        }
        return tag;
    }

    public void readTagPair(String infileS) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
            int tagNum = dis.readInt();
            this.tagLengthInLong = dis.readInt();
            this.tags = new long[tagLengthInLong][tagNum];
            this.tagLength = new byte[tagNum];
            this.order = new int[tagNum];
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = dis.readLong();
                }
                tagLength[i] = dis.readByte();
                order[i] = dis.readInt();
            }
            System.out.println(tagNum / 2 + " TagPairs are read");
        } catch (Exception e) {
            System.out.println("Error occurred while reading " + infileS + " " + e.toString());
        }
    }

    public void writeTagPair(String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tags[0].length);
            dos.writeInt(tagLengthInLong);
            for (int i = 0; i < tags[0].length; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    dos.writeLong(tags[j][i]);
                }
                dos.writeByte(tagLength[i]);
                dos.writeInt(order[i]);
            }
            dos.flush();
            dos.close();
            System.out.println(tags[0].length / 2 + " TagPairs are written");
        } catch (Exception e) {
            System.out.println("Error occurred while writing " + outfileS + " " + e.toString());
        }
    }

    public void sortBySeq() {
        GenericSorting.quickSort(0, tags[0].length, compSeq, swapper);
        System.out.println("TagPair is sorted by sequence");
    }

    public void sortByOrder() {
        GenericSorting.quickSort(0, tags[0].length, compOrder, swapper);
        System.out.println("TagPair is sorted by pair order");
    }
    Swapper swapper = new Swapper() {

        public void swap(int a, int b) {
            long t1, t2;
            for (int i = 0; i < tagLengthInLong; i++) {
                t1 = tags[i][a];
                tags[i][a] = tags[i][b];
                tags[i][b] = t1;
            }
            byte t3;
            t3 = tagLength[a];
            tagLength[a] = tagLength[b];
            tagLength[b] = t3;
            int t4;
            t4 = order[a];
            order[a] = order[b];
            order[b] = t4;
        }
    };
    IntComparator compSeq = new IntComparator() {

        public int compare(int a, int b) {
            for (int i = 0; i < tagLengthInLong; i++) {
                if (tags[i][a] < tags[i][b]) {
                    return -1;
                }
                if (tags[i][a] > tags[i][b]) {
                    return 1;
                }
            }
            return 0;
        }
    };
    IntComparator compOrder = new IntComparator() {

        public int compare(int a, int b) {
            return order[a] < order[b] ? -1 : 1;
        }
    };
}
