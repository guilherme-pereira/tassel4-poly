/*
 * Clusters
 */
package net.maizegenetics.gbs.pipeline;

import java.io.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import net.maizegenetics.gbs.homology.TagMatchFinder;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.gbs.util.ReadsByTaxa;

/**
 *
 * @author fl262
 */
public class Clusters {

    cluster[] cls;

    public Clusters(ReadsByTaxa rbt) {
        getClusters(rbt);
    }

    public Clusters(TagsByTaxa tbt) {
        getClusters(tbt);
    }

    public Clusters(String infileS, boolean binary) {
        readCluster(infileS, binary);
    }

    public void getClusters(ReadsByTaxa rbt) {
        PolymorphismFinder pf = new PolymorphismFinder(rbt);
        ArrayList<cluster> clList = new ArrayList();
        for (int i = 0; i < rbt.haplotypeNum; i++) {
            long[] queryLongSeq = new long[2];
            queryLongSeq[0] = rbt.haplotype[0][i];
            queryLongSeq[1] = rbt.haplotype[1][i];
            ArrayList<Integer> hitIndex = pf.findOneMismatch(queryLongSeq);
            if (hitIndex.isEmpty()) {
                continue;
            }
            Integer[] hitIndexArray = hitIndex.toArray(new Integer[hitIndex.size()]);
            for (int j = 0; j < hitIndexArray.length; j++) {
                cluster cl = new cluster(i, hitIndexArray[j], hitIndexArray.length + 1, true);
                clList.add(cl);
            }
        }
        cls = clList.toArray(new cluster[clList.size()]);
    }

    // This constructor works for tbt file
    public void getClusters(TagsByTaxa tbt) {
        ArrayList<cluster> clList = new ArrayList();
        TagMatchFinder tmf = new TagMatchFinder(tbt);
        for (int i = 0; i < tbt.getTagCount(); i++) {
            long[] qTag = tbt.getTag(i);
            TreeMap<Integer, Integer> hitDiv = tmf.findMatchesWithIntLengthWords(qTag, 1, false);
            for (Entry<Integer, Integer> each : hitDiv.entrySet()) {
                if (each.getValue() > 0) {
                    clList.add(new cluster(i, each.getKey(), hitDiv.size(), true));
                }
            }
        }
        cls = clList.toArray(new cluster[clList.size()]);
        //Arrays.sort(cls);
    }

    public void networkFilter() {
        TreeSet<cluster> clSet = new TreeSet();
        TreeSet<Integer> threeSet = new TreeSet();
        ArrayList<cluster> twoClusterList = new ArrayList();
        cluster[] twoCluster;
        for (int i = 0; i < cls.length; i++) {
            cls[i].switchQueryAndHit();
        }
        Arrays.sort(cls);
        for (int i = 0; i < cls.length; i++) {
            if (cls[i].cSize > 2) {
                threeSet.add(cls[i].queryIndex);
                threeSet.add(cls[i].hitIndex);
            } else {
                clSet.add(cls[i]);
            }
        }
        cluster[] tempTwoCluster = clSet.toArray(new cluster[clSet.size()]);
        Integer[] threeArray = threeSet.toArray(new Integer[threeSet.size()]);
        Arrays.sort(threeArray);
        for (int i = 0; i < tempTwoCluster.length; i++) {
            int hit = Arrays.binarySearch(threeArray, tempTwoCluster[i].queryIndex);
            if (hit > -1) {
                continue;
            }
            hit = Arrays.binarySearch(threeArray, tempTwoCluster[i].hitIndex);
            if (hit > -1) {
                continue;
            }
            twoClusterList.add(tempTwoCluster[i]);
        }
        twoCluster = twoClusterList.toArray(new cluster[twoClusterList.size()]);
        cls = twoCluster;
        Arrays.sort(cls);
        twoCluster = null;
        tempTwoCluster = null;
        threeArray = null;
        clSet = null;
        threeSet = null;
        twoClusterList = null;
    }

    public void repeatLibraryFilter(ReadsByTaxa rbt, String TagSeqS, String repeatlibS) {
        String[] tagSeqs = new String[this.getSnpCount()];
        int[] clsIndex = new int[tagSeqs.length];
        int count = 0;
        for (int i = 0; i < cls.length; i++) {
            if (!cls[i].ifSnp) {
                continue;
            }
            long[] longSeq = new long[2];
            longSeq[0] = rbt.haplotype[0][cls[i].queryIndex];
            longSeq[1] = rbt.haplotype[1][cls[i].queryIndex];
            tagSeqs[count] = BaseEncoder.getSequenceFromLong(longSeq);
            clsIndex[count] = i;
            count++;
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(TagSeqS), 65536);
            for (int i = 0; i < tagSeqs.length; i++) {
                bw.write(">" + String.valueOf(clsIndex[i]));
                bw.newLine();
                bw.write(tagSeqs[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error occurred while writing " + TagSeqS);
        }
        String blastFileS = TagSeqS.replace("fasta", "blast");
        String cmd = "blastn -query " + TagSeqS + " -db " + repeatlibS + " -out " + blastFileS + " -evalue 1e-10 -outfmt 6";
        System.out.println(cmd);
        try {
            Runtime.getRuntime().exec(cmd).waitFor();
        } catch (Exception e) {
            System.out.println("Error occurred while doing blast " + e.toString());
        }
        TreeSet<Integer> repeatIndex = new TreeSet();
        try {
            BufferedReader br = new BufferedReader(new FileReader(blastFileS), 65536);
            String temp;
            while ((temp = br.readLine()) != null) {
                Integer index = Integer.valueOf(temp.split("\t")[0]);
                repeatIndex.add(index);
            }
        } catch (Exception e) {
            System.out.println("Error occurred while reading" + blastFileS);
        }
        Integer[] repeatIndexArray = repeatIndex.toArray(new Integer[repeatIndex.size()]);
        for (int i = 0; i < repeatIndexArray.length; i++) {
            cls[repeatIndexArray[i]].ifSnp = false;
        }
    }

    public void repeatFilter(ReadsByTaxa rbt, double ratio) {
        int[] tagCount = new int[2 * this.getSnpCount()];
        int[] clusterIndex = new int[2 * this.getSnpCount()];
        for (int i = 0; i < cls.length; i++) {
            if (!cls[i].ifSnp) {
                continue;
            }
            int totalA = 0, totalC = 0;
            for (int j = 0; j < rbt.taxaNum; j++) {
                totalA += rbt.hapDist[cls[i].queryIndex][j];
                totalC += rbt.hapDist[cls[i].hitIndex][j];
            }
            tagCount[2 * i] = totalA;
            clusterIndex[2 * i] = i;
            tagCount[2 * i + 1] = totalC;
            clusterIndex[2 * i + 1] = i;
        }
        for (int i = 0; i < tagCount.length - 1; i++) {
            for (int j = i + 1; j < tagCount.length; j++) {
                if (tagCount[i] > tagCount[j]) {
                    int mid = tagCount[i];
                    tagCount[i] = tagCount[j];
                    tagCount[j] = mid;
                    mid = clusterIndex[i];
                    clusterIndex[i] = clusterIndex[j];
                    clusterIndex[j] = mid;
                }
            }
        }
        int beginIndex = (int) Math.floor(clusterIndex.length * (1 - ratio));
        for (int i = beginIndex; i < clusterIndex.length; i++) {
            cls[clusterIndex[i]].ifSnp = false;
        }
    }

    public int getSnpCount() {
        int snpCount = 0;
        for (int i = 0; i < cls.length; i++) {
            if (cls[i].ifSnp) {
                snpCount++;
            }
        }
        return snpCount;
    }

    public void alleleFrequencyFileter(ReadsByTaxa rbt, double minBorderMaf, double maxBorderMaf) {
        for (int i = 0; i < cls.length; i++) {
            if (!cls[i].ifSnp) {
                continue;
            }
            int totalA = 0, totalC = 0;
            for (int j = 0; j < rbt.taxaNum; j++) {
                totalA += rbt.hapDist[cls[i].queryIndex][j];
                totalC += rbt.hapDist[cls[i].hitIndex][j];
            }
            int min = totalA;
            if (totalC < totalA) {
                min = totalC;
            }
            double maf = (double) min / (double) (totalA + totalC);
            if (maf < minBorderMaf || maf > maxBorderMaf) {
                cls[i].ifSnp = false;
            }
        }
    }

    public void heteozygoteFilter(ReadsByTaxa rbt) {
        for (int i = 0; i < cls.length; i++) {
            int heteoA = 0, heteoC = 0, totalA = 0, totalC = 0;
            boolean flag = false;
            if (cls[i].ifSnp) {
                for (int j = 0; j < rbt.taxaNum; j++) {
                    if (rbt.hapDist[cls[i].queryIndex][j] > 0) {
                        totalA += rbt.hapDist[cls[i].queryIndex][j];
                        if (rbt.hapDist[cls[i].hitIndex][j] > 0) {
                            totalC += rbt.hapDist[cls[i].hitIndex][j];
                            heteoA += rbt.hapDist[cls[i].queryIndex][j];
                            heteoC += rbt.hapDist[cls[i].hitIndex][j];
                            boolean ifSig = chiSquareEvenDf1(rbt.hapDist[cls[i].queryIndex][j], rbt.hapDist[cls[i].hitIndex][j], 3.841);
                            if (ifSig) {
                                cls[i].ifSnp = false;
                                flag = true;
                                break;
                            }
                        }
                    } else {
                        if (rbt.hapDist[cls[i].hitIndex][j] > 0) {
                            totalC += rbt.hapDist[cls[i].hitIndex][j];
                        }
                    }
                }
            } else {
                continue;
            }
            if (flag) {
                continue;
            }


            if (chiSquareEvenDf1(heteoA, heteoC, 3.841)) {
                cls[i].ifSnp = false;
                continue;
            }
            /*
            if (!chiSquareDf1(heteoA, heteoC, totalA, totalC, 0)) {
            cls[i].ifSnp = false;
            }
             * 
             */

        }
    }

    public boolean chiSquareDf1(int ob1, int ob2, int ex1, int ex2, double chiValueCutoff) {
        double o1 = (double) ob1 / (double) (ob1 + ob2);
        double o2 = (double) ob2 / (double) (ob1 + ob2);
        double e1 = (double) ex1 / (double) (ex1 + ex2);
        double e2 = (double) ex2 / (double) (ex1 + ex2);
        double chi = Math.pow(o1 - e1, 2) / e1 + Math.pow(o2 - e2, 2) / e2;
        if (chi > chiValueCutoff) {
            return true;
        }
        return false;
    }

    public boolean chiSquareEvenDf1(int ob1, int ob2, double chiValueCutoff) {
        double o1 = (double) ob1;
        double o2 = (double) ob2;
        double ex = (o1 + o2) / 2;
        double chi = Math.pow(o1 - ex, 2) / ex + Math.pow(o2 - ex, 2) / ex;
        if (chi > chiValueCutoff) {
            return true;
        }
        return false;
    }

    public void readCluster(String infileS, boolean binary) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
            cls = new cluster[dis.readInt()];
            int[] temp = new int[3];
            for (int i = 0; i < cls.length; i++) {
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = dis.readInt();
                }
                cls[i] = new cluster(temp[0], temp[1], temp[2], dis.readBoolean());
            }
        } catch (Exception e) {
            System.out.println("erroe in reading" + infileS + e.toString());
        }
    }

    public void hapDetail(ReadsByTaxa rbt, String outfileS, float coverRate) {
        ArrayList<String> al = new ArrayList();
        String[] hapMap;
        for (int i = 0; i < cls.length; i++) {
            if (cls[i].ifSnp) {
                String temp;
                temp = getHapDetail(rbt, i, cls[i].queryIndex, cls[i].hitIndex, coverRate);
                if (temp != null) {
                    al.add(temp);
                }
            }
        }
        hapMap = al.toArray(new String[al.size()]);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
            bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
            for (int i = 0; i < rbt.getTaxaCount() - 1; i++) {
                bw.write(rbt.getTaxaNames()[i] + "\t");
            }
            bw.write(rbt.getTaxaNames()[rbt.getTaxaCount() - 1]);
            bw.newLine();
            int count = 0;
            for (int i = 0; i < hapMap.length; i++) {
                bw.write(hapMap[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println(e.toString());
        }
    }

    private String getHapDetail(ReadsByTaxa rbt, int ID, int queryIndex, int hitIndex, float coverRate) {
        StringBuilder sb = new StringBuilder();
        sb.append(queryIndex).append("\tA/C\t").append(1).append("\t").append(hitIndex).append("\t+\tNA\tSWGDiv\tGBS\tSWGV1\tSWGPop\tQC+\t");
        int countN = 0;
        int totalA = 0, totalC = 0;
        int hetQuery = 0, hetHit = 0;
        for (int i = 0; i < rbt.taxaNum; i++) {
            if (rbt.hapDist[queryIndex][i] > 0) {
                totalA += rbt.hapDist[queryIndex][i];
                if (rbt.hapDist[hitIndex][i] > 0) {
                    totalC += rbt.hapDist[hitIndex][i];
                    sb.append(String.valueOf(rbt.hapDist[queryIndex][i])).append("|").append(String.valueOf(rbt.hapDist[hitIndex][i])).append("\t");
                    hetQuery += rbt.hapDist[queryIndex][i];
                    hetHit += rbt.hapDist[hitIndex][i];
                } else {
                    sb.append(String.valueOf(rbt.hapDist[queryIndex][i])).append("|").append("\t");
                }
            } else {
                if (rbt.hapDist[hitIndex][i] > 0) {
                    totalC += rbt.hapDist[hitIndex][i];
                    sb.append("|");
                    sb.append(String.valueOf(rbt.hapDist[hitIndex][i])).append("\t");
                } else {
                    sb.append("N").append("\t");
                    countN++;
                }
            }
        }

        if (1 - (float) countN / (float) rbt.taxaNum < coverRate) {
            return null;
        }
        sb.append(String.valueOf(hetQuery)).append("\t");
        sb.append(String.valueOf(hetHit)).append("\t");
        sb.append(String.valueOf(totalA)).append("\t");
        sb.append(String.valueOf(totalC));
        return sb.toString();
    }

    public void writeHapMap(ReadsByTaxa rbt, String outfileS, float coverRate) {
        String[] hapMap;
        ArrayList<String> al = new ArrayList();
        for (int i = 0; i < cls.length; i++) {
            String temp;
            if (cls[i].ifSnp) {
                temp = getHapMapRecord(rbt, i + 1, cls[i].queryIndex, cls[i].hitIndex, coverRate);
                if (temp != null) {
                    al.add(temp);
                }
            }
        }
        hapMap = al.toArray(new String[al.size()]);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
            bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
            for (int i = 0; i < rbt.getTaxaCount() - 1; i++) {
                bw.write(rbt.getTaxaNames()[i] + "\t");
            }
            bw.write(rbt.getTaxaNames()[rbt.getTaxaCount() - 1]);
            bw.newLine();
            for (int i = 0; i < hapMap.length; i++) {
                bw.write(hapMap[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println(e.toString());
        }
    }

    private String getHapMapRecord(ReadsByTaxa rbt, int ID, int queryIndex, int hitIndex, float coverRate) {
        StringBuilder sb = new StringBuilder();
        sb.append(queryIndex).append("|").append(hitIndex).append("\tA/C\t").append(1).append("\t").append(ID).append("\t+\tNA\tSWGDiv\tGBS\tSWGV1\tSWGPop\tQC+\t");
        int countN = 0;
        int totalA = 0, totalC = 0;
        for (int i = 0; i < rbt.taxaNum; i++) {
            if (rbt.hapDist[queryIndex][i] > 0) {
                totalA += rbt.hapDist[queryIndex][i];
                if (rbt.hapDist[hitIndex][i] > 0) {
                    totalC += rbt.hapDist[hitIndex][i];
                    sb.append("R").append("\t");
                } else {
                    sb.append("A").append("\t");
                }
            } else {
                if (rbt.hapDist[hitIndex][i] > 0) {
                    totalC += rbt.hapDist[hitIndex][i];
                    sb.append("C").append("\t");
                } else {
                    sb.append("N").append("\t");
                    countN++;
                }
            }
        }

        if (1 - (float) countN / (float) rbt.taxaNum < coverRate) {
            return null;
        }
        sb.deleteCharAt(sb.length() - 1);
        return sb.toString();
    }

    public void writeCluster(String outfileS, boolean binary) {
        int snpCount = getSnpCount();
        if (binary) {
            try {
                DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
                dos.writeInt(snpCount);
                for (int i = 0; i < cls.length; i++) {
                    cls[i].writeBinary(dos);
                }
                dos.flush();
                dos.close();
            } catch (Exception e) {
                System.out.println(e.toString());
            }
        } else {
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
                bw.write(Integer.toString(snpCount));
                bw.newLine();
                bw.write("queryIndex\thitIndex\tdiv\tclusterSize");
                bw.newLine();
                for (int i = 0; i < cls.length; i++) {
                    cls[i].writeTxt(bw);
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                System.out.println(e.toString());
            }
        }
    }

    public void writeFastA(ReadsByTaxa rbt, String outFastAS) {
        try {
            int count = 1;
            BufferedWriter bw = new BufferedWriter(new FileWriter(outFastAS), 65536);
            for (int i = 0; i < cls.length; i++) {
                if (!cls[i].ifSnp) {
                    continue;
                }
                long longSeq[] = new long[2];
                longSeq[0] = rbt.haplotype[0][cls[i].queryIndex];
                longSeq[1] = rbt.haplotype[1][cls[i].queryIndex];
                bw.write(">" + String.valueOf(count));
                bw.newLine();
                bw.write("G" + BaseEncoder.getSequenceFromLong(longSeq));
                bw.newLine();
                count++;
                longSeq[0] = rbt.haplotype[0][cls[i].hitIndex];
                longSeq[1] = rbt.haplotype[1][cls[i].hitIndex];
                bw.write(">" + String.valueOf(count));
                bw.newLine();
                bw.write("G" + BaseEncoder.getSequenceFromLong(longSeq));
                bw.newLine();
                count++;
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error occurred while writing" + outFastAS);
        }
    }

    public void prin() {
        for (int i = 0; i < 100; i++) {
            cls[i].screenPri();
        }
        System.out.println(cls.length);
    }
}

class cluster implements Comparable<cluster> {

    int queryIndex;
    int hitIndex;
    int cSize;
    boolean ifSnp;

    public cluster(int queryIndex, int hitIndex, int cSize, boolean ifSnp) {
        this.queryIndex = queryIndex;
        this.hitIndex = hitIndex;
        this.cSize = cSize;
        this.ifSnp = ifSnp;
    }

    public void setIfSnp(boolean ifSnp) {
        this.ifSnp = ifSnp;
    }

    public void writeBinary(DataOutputStream dos) {
        try {
            if (ifSnp) {
                dos.writeInt(queryIndex);
                dos.writeInt(hitIndex);
                dos.writeInt(cSize);
                dos.writeBoolean(ifSnp);
            }
        } catch (Exception e) {
            System.out.println(e.toString());
        }
    }

    public void writeTxt(BufferedWriter bw) {
        try {
            if (ifSnp) {
                bw.write(queryIndex + "\t" + hitIndex + "\t" + cSize + "\t" + String.valueOf(ifSnp));
                bw.newLine();
            }
        } catch (Exception e) {
            System.out.println(e.toString());
        }
    }

    public void screenPri() {
        System.out.println(queryIndex + "\t" + hitIndex + "\t" + cSize);
    }

    public void switchQueryAndHit() {
        if (queryIndex > hitIndex) {
            int mid = queryIndex;
            queryIndex = hitIndex;
            hitIndex = mid;
        }
    }

    public int compareTo(cluster o) {
        if (queryIndex < o.queryIndex) {
            return -1;
        } else if (queryIndex > o.queryIndex) {
            return 1;
        } else {
            if (hitIndex < o.hitIndex) {
                return -1;
            } else if (hitIndex > o.hitIndex) {
                return 1;
            } else {
                if (cSize < o.cSize) {
                    return -1;
                } else if (cSize > o.cSize) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }
    }
}
