/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pd;
import java.io.File;
import java.util.HashSet;

/**
* User: dkroon
*/
public class FolderParser {

    private File directory;
    private int traitPosition = 0;
    private int chrPosition = 5;
    private String filenameToken = "\\.";

    public FolderParser(File directoryIn){
        this.directory = directoryIn;

    }

    public FolderParser(File directoryIn, int traitPositionIn, int chrPositionIn, String filenameTokenIn){
        if(traitPosition >= 0) traitPosition = traitPositionIn;
        if(chrPosition >= 0) chrPosition = chrPositionIn;
        if(filenameToken!= null) filenameToken = filenameTokenIn;
    }

    public String[] getAllTraits(){

        File[] files = getFiles();

        HashSet<String> traitSet = new HashSet<String>();
        for(int i=0; i<files.length; i++){
            traitSet.add(getFileFeature(files[i],filenameToken, traitPosition));
        }
        String[] trait = new String[traitSet.size()];
        traitSet.toArray(trait);
        return trait;
    }

    public File getFile(String chromosome, String trait){
        File[] chrFiles = getAllChromosomeFiles(chromosome);
        for (int i = 0; i < chrFiles.length ; i++) {
            if(chrFiles[i].toString().contains(trait)) return chrFiles[i];
        }
        return null;
    }

    private File[] getAllChromosomeFiles(String chromosome){

        HashSet<File> chrFile = new HashSet<File>();

        File[] files = getFiles();
        for(int i = 0; i < files.length; i++){
            if(files[i].toString().contains(chromosome)) chrFile.add(files[i]);
        }

        File[] chrFiles = new File[chrFile.size()];
        chrFile.toArray(chrFiles);
        return chrFiles;
    }


    private String getFileFeature(File aFile, String token, int field){
        String fileName = aFile.getName();
        String[] part = fileName.split(token);
        return part[field];
    }

    private File[] getFiles(){
        File[] directoryFile = directory.listFiles();
        if (directoryFile.length < 1) {
            String output = "Directory contains no files.";
// JOptionPane.showMessageDialog(c, "Directory contains no files: \n" + directory.getPath());
            return null;
        }
        return directoryFile;
    }



    public static void main(String[] args) {

        String testDir = "C:\\Users\\dkroon\\Documents\\TASSEL\\code\\GWAS\\src\\supporting\\gwas_results";
        File fileDir = new File(testDir);
        FolderParser fp = new FolderParser(fileDir);

        File[] file = fp.getAllChromosomeFiles("chr8");

        System.out.println("numberOfChrFiles = " + file.length);
        for (int i = 0; i < file.length ; i++) {
            System.out.println("file = " + file[i]);
        }

        String[] trait = fp.getAllTraits();

        System.out.println("numberOfTraitsAcrossAllFiles = " + trait.length);
        for (int i = 0; i < trait.length; i++) {
            System.out.println("trait = " + trait[i]);
        }
    }
}