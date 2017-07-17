package net.maizegenetics.gbs.maps;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Update a physical map with genetic mapping information.
 *
 * The genetic mapping results would normally be derived from a run of
 * TagCallerAgainstAnchorMT.
 *
 * @author edbuckler
 */
public class UpdateTOPMWithGeneticMapping {

    TagsOnPhysicalMap theTOPM = null;
    double minRatioOfBestToNext = 2.0;
    int maxDistBLASTLD = 10000000;

    public UpdateTOPMWithGeneticMapping(String sTOPM, String sMappings, String sOutFile) {
        theTOPM = new TagsOnPhysicalMap(sTOPM, sTOPM.contains(".bin"));
        try {
            BufferedReader br = new BufferedReader(new FileReader(sMappings));
            String temp;
            int chrLD, positionLD, countBetterSites, tagCount;
            double ratioOfBestToNext, bestPValue;
            int count = 0, revisedGeneticPositions = 0, newGeneticPositions = 0, agreeWithBLAST = 0;
            int hasBLASTPosition = 0, geneticallyMapsSig1e7 = 0, tagsWGeneticOrBLASTSig = 0;
            long readCountTotal = 0, readCountMap = 0, readCountNotMap = 0, readWGeneticOrBLASTSig = 0;
            while ((temp = br.readLine()) != null) {
                String[] s = temp.split("\\s");
                //               System.out.println(Arrays.toString(s));
                if (s.length != 17) {
                    continue;
                }
                try {
                    chrLD = Integer.parseInt(s[5]);
                    positionLD = Integer.parseInt(s[7]);
                    countBetterSites = Integer.parseInt(s[14]);
                    ratioOfBestToNext = Double.parseDouble(s[12]);
                    bestPValue = Double.parseDouble(s[8]);
                    tagCount = Integer.parseInt(s[10]);
                    //           System.out.printf("%d %d %d %g %g %n",chrLD,positionLD,countBetterSites,ratioOfBestToNext,bestPValue);
                } catch (NumberFormatException nef) {
                    System.err.println("Error in Reading Lines:" + temp);
                    continue;
                }
                count++;
                long[] tag = BaseEncoder.getLongArrayFromSeq(s[0]);
                int tagIndex = theTOPM.getTagIndex(tag);
                if (tagIndex < 0) {
                    continue;
                }
                readCountTotal += tagCount;
                if (theTOPM.getChromosome(tagIndex) > -1) {
                    hasBLASTPosition++;
                }
                if (bestPValue < 1e-7) {
                    geneticallyMapsSig1e7++;
                }
                if ((theTOPM.getChromosome(tagIndex) > -1) || (bestPValue < 1e-7) || ((ratioOfBestToNext > minRatioOfBestToNext) && (countBetterSites > 2))) {
                    tagsWGeneticOrBLASTSig++;
                    readWGeneticOrBLASTSig += tagCount;
                }
                int distance = Math.abs(positionLD - theTOPM.getStartPosition(tagIndex));
                if ((chrLD == theTOPM.getChromosome(tagIndex)) && (distance < maxDistBLASTLD)) {
                    theTOPM.setMapP(tagIndex, bestPValue);
                    agreeWithBLAST++;
                    readCountMap += tagCount;
                } else if ((chrLD != theTOPM.getChromosome(tagIndex)) && (ratioOfBestToNext > minRatioOfBestToNext) && (countBetterSites > 2)) {
                    if (theTOPM.getChromosome(tagIndex) < 0) {
                        newGeneticPositions++;
                    } else {
                        revisedGeneticPositions++;
                    }
                    updateTOPM(tag, chrLD, positionLD, positionLD,
                            TagsOnPhysicalMap.BYTE_MISSING, TagsOnPhysicalMap.BYTE_MISSING, bestPValue);
                    readCountMap += tagCount;
                } else {
                    //                   System.out.println(temp);
                    readCountNotMap += tagCount;
                }
                if (count % 10000 == 0) {
                    System.out.printf("Count: %d  hasBLAST: %d GeneticLT1e-7: %d hasBLASTorGeneticSig: %d NewPos: %d  RevisedPos: %d ConfirmedPos: %d"
                            + " TotalReads: %d ReadWithBLASTorGenetics: %d TotalReadMapCount: %d TotalReadNoMapCount: %d %n",
                            count, hasBLASTPosition, geneticallyMapsSig1e7, tagsWGeneticOrBLASTSig, newGeneticPositions, revisedGeneticPositions, agreeWithBLAST,
                            readCountTotal, readWGeneticOrBLASTSig, readCountMap, readCountNotMap);
                }
//                if(count%100==0) {
//                    System.out.println(temp);
//                    System.out.println(theTOPM.printRow(tagIndex));
//                }
            }
//                           if(theTOPM==null) continue;
            System.out.printf("Count: %d  hasBLAST: %d GeneticLT1e-7: %d hasBLASTorGeneticSig: %d NewPos: %d  RevisedPos: %d ConfirmedPos: %d"
                    + " TotalReads: %d ReadWithBLASTorGenetics: %d TotalMapCount: %d NoMapCount: %d %n",
                    count, hasBLASTPosition, geneticallyMapsSig1e7, tagsWGeneticOrBLASTSig, newGeneticPositions, revisedGeneticPositions, agreeWithBLAST,
                    readCountTotal, readWGeneticOrBLASTSig, readCountMap, readCountNotMap);
            theTOPM.printRows(100000, true, 10);
//            theTOPM.writeTextFile(new File(sOutFile));
//            theTOPM.writeBinaryFile(new File(sOutFile));
//                if(i%10000==0) theTOPM.writeBinaryFile(new File("/Users/edbuckler/SolexaAnal/GBS/test/14FCGBS_gen110404.tg.ndup.bin"));
//
        } catch (IOException e) {
            System.err.println("Error in UpdateTOPMWithGeneticMapping");

        }

    }

    private void updateTOPM(long[] tag, int chr, int minPosition, int maxPosition, byte strand, byte divergence, double mapP) {
        int index = theTOPM.getTagIndex(tag);
        if (index > -1) {
            theTOPM.setChromoPosition(index, chr, strand, minPosition, maxPosition);
            theTOPM.setDivergence(index, divergence);
            theTOPM.setMapP(index, mapP);
        }
    }

    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println("Input format -topm MapFile -map Mapping -o OutMapFile");
//            String[] s={"-topm","/Users/edbuckler/SolexaAnal/GBS/reftags/14FCGBS.tg.ndup.bin",
//            "-map","/Users/edbuckler/SolexaAnal/GBS/mapping/allfusionMapping110408.txt",
//            "-o","/Users/edbuckler/SolexaAnal/GBS/test/allfusionMapping110408x.topm.txt"};
            String[] s = {"-topm", "/Volumes/LaCie/mergedNam282Ames_072011_mappedOnly.topm.bin",
                "-map", "/Users/edbuckler/SolexaAnal/GBS/build111217/mapping/Zea20111217_mapping_chr10.txt",
                "-o", "/Users/edbuckler/SolexaAnal/GBS/build111217/test/C10mergedNam282Ames_072011_mappedOnly.topm.bin"};
            args = s;
        }
        if (args.length != 6) {
            System.err.println("Input format -topm MapFile -map Mapping -o OutMapFile");
            return;
        }
        UpdateTOPMWithGeneticMapping theUTGM = new UpdateTOPMWithGeneticMapping(args[1], args[3], args[5]);
    }
}
