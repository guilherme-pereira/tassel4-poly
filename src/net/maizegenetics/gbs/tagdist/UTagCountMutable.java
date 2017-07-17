/*
 * UTagCountMutable
 */
package net.maizegenetics.gbs.tagdist;

import java.util.ArrayList;

/**
 *
 * @author Fei Lu
 */
public class UTagCountMutable extends TagCountMutable {

    ArrayList<Long>[] tagsList;
    ArrayList<Byte> tagLengthList;
    ArrayList<Integer> readCountList;

    public UTagCountMutable(int tagLengthInLong, int maxSize) {
        super(tagLengthInLong, 0);
        tagsList = new ArrayList[tagLengthInLong];
        for (int i = 0; i < tagLengthInLong; i++) {
            tagsList[i] = new ArrayList();
        }
        tagLengthList = new ArrayList();
        readCountList = new ArrayList();
    }

    @Override
    public void addReadCount(long[] tag, int length, int count) {
        for (int i = 0; i < tag.length; i++) {
            tagsList[i].add(tag[i]);
        }
        tagLengthList.add((byte) length);
        readCountList.add(count);
        currentRows++;
    }

    public void toArray() {
        int tagNum = tagLengthList.size();
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        readCount = new int[tagNum];
        for (int i = 0; i < tagNum; i++) {
            for (int j = 0; j < tagsList.length; j++) {
                tags[j][i] = tagsList[j].get(i);
            }
            tagLength[i] = tagLengthList.get(i);
            readCount[i] = readCountList.get(i);
        }
    }
}
