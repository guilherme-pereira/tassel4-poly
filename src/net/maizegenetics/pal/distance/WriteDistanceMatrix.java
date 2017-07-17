package net.maizegenetics.pal.distance;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;


/**
 * @author Terry
 */
public class WriteDistanceMatrix {

    private WriteDistanceMatrix() {
    }

    public static void saveDelimitedDistanceMatrix(DistanceMatrix matrix, String saveFile) {
        saveDelimitedDistanceMatrix(matrix, new File(saveFile));
    }

    public static void saveDelimitedDistanceMatrix(DistanceMatrix matrix, File saveFile) {

        if (saveFile == null) {
            return;
        }
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {

            fw = new FileWriter(saveFile);
            bw = new BufferedWriter(fw);

            bw.write(String.valueOf(matrix.getRowCount()));
            bw.write("\n");

            for (int r = 0, n = matrix.getRowCount(); r < n; r++) {
                Object[] theRow = matrix.getRow(r);
                for (int i = 0; i < theRow.length; i++) {
                    if (i != 0) {
                        bw.write("\t");
                    }
                    bw.write(theRow[i].toString());
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            System.out.println("WriteDistanceMatrix: saveDelimitedDistanceMatrix: problem writing file: " + e.getMessage());
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }
}
