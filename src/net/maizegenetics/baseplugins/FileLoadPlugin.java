/*
 * FileLoadPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.ReadDistanceMatrix;
import net.maizegenetics.pal.report.Report;
import net.maizegenetics.pal.report.TableReportUtils;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import net.maizegenetics.prefs.TasselPrefs;

import net.maizegenetics.gui.DialogUtils;

import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import javax.swing.*;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.io.*;

import java.net.URL;

import java.util.List;
import java.util.ArrayList;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

/**
 *
 * @author Ed Buckler
 */
public class FileLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FileLoadPlugin.class);
    private String[] myOpenFiles = null;
    private TasselFileType myFileType = TasselFileType.Unknown;
    private PlinkLoadPlugin myPlinkLoadPlugin = null;
    private FlapjackLoadPlugin myFlapjackLoadPlugin = null;
    private boolean myIsSBit = true;
    private JFileChooser myOpenFileChooser = new JFileChooser(TasselPrefs.getOpenDir());

    public enum TasselFileType {

        SqrMatrix, Annotated, Sequence, Polymorphism, Numerical, Unknown, Fasta,
        Hapmap, Plink, Phenotype, Flapjack, Phylip_Seq, Phylip_Inter, GeneticMap, Table,
        Serial, HapmapDiploid, Text, HDF5, VCF, ByteHDF5
    };
    public static final String FILE_EXT_HAPMAP = ".hmp.txt";
    public static final String FILE_EXT_HAPMAP_GZ = ".hmp.txt.gz";
    public static final String FILE_EXT_PLINK_MAP = ".plk.map";
    public static final String FILE_EXT_PLINK_PED = ".plk.ped";
    public static final String FILE_EXT_FLAPJACK_MAP = ".flpjk.map";
    public static final String FILE_EXT_FLAPJACK_GENO = ".flpjk.geno";
    public static final String FILE_EXT_SERIAL_GZ = ".serial.gz";
    public static final String FILE_EXT_HDF5 = ".hmp.h5";
    public static final String FILE_EXT_VCF = ".vcf";

    /**
     * Creates a new instance of FileLoadPlugin
     */
    public FileLoadPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public FileLoadPlugin(Frame parentFrame, boolean isInteractive, PlinkLoadPlugin plinkLoadPlugin, FlapjackLoadPlugin flapjackLoadPlugin) {
        super(parentFrame, isInteractive);
        myPlinkLoadPlugin = plinkLoadPlugin;
        myFlapjackLoadPlugin = flapjackLoadPlugin;
    }

    public DataSet performFunction(DataSet input) {

        try {

            if (isInteractive()) {
                FileLoadPluginDialog theDialog = new FileLoadPluginDialog();
                theDialog.setLocationRelativeTo(getParentFrame());
                theDialog.setVisible(true);
                if (theDialog.isCancel()) {
                    return null;
                }
                myFileType = theDialog.getTasselFileType();

                if (myFileType == TasselFileType.Plink) {
                    return myPlinkLoadPlugin.performFunction(null);
                }

                if (myFileType == TasselFileType.Flapjack) {
                    return myFlapjackLoadPlugin.performFunction(null);
                }

                setOpenFiles(getOpenFilesByChooser());
                theDialog.dispose();
            }

            if ((myOpenFiles == null) || (myOpenFiles.length == 0)) {
                return null;
            }

            List result = new ArrayList();
            ArrayList<String> alreadyLoaded = new ArrayList();
            for (int i = 0; i < myOpenFiles.length; i++) {

                if (alreadyLoaded.contains(myOpenFiles[i])) {
                    continue;
                }

                DataSet tds = null;
                try {

                    if (myFileType == TasselFileType.Unknown) {
                        if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP) || myOpenFiles[i].endsWith(FILE_EXT_HAPMAP_GZ)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Hapmap);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_PED)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                            String theMapFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_PED, FILE_EXT_PLINK_MAP);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(theMapFile);
                            myPlinkLoadPlugin.loadFile(myOpenFiles[i], theMapFile, null);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_MAP)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                            String thePedFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_MAP, FILE_EXT_PLINK_PED);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(thePedFile);
                            myPlinkLoadPlugin.loadFile(thePedFile, myOpenFiles[i], null);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_FLAPJACK_GENO)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Flapjack);
                            String theMapFile = myOpenFiles[i].replaceFirst(FILE_EXT_FLAPJACK_GENO, FILE_EXT_FLAPJACK_MAP);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(theMapFile);
                            myFlapjackLoadPlugin.loadFile(myOpenFiles[i], theMapFile, null);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_FLAPJACK_MAP)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Flapjack);
                            String theGenoFile = myOpenFiles[i].replaceFirst(FILE_EXT_FLAPJACK_MAP, FILE_EXT_FLAPJACK_GENO);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(theGenoFile);
                            myFlapjackLoadPlugin.loadFile(theGenoFile, myOpenFiles[i], null);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_SERIAL_GZ)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Serial);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Serial);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_HDF5)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.HDF5);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.HDF5);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_VCF) || myOpenFiles[i].endsWith(FILE_EXT_VCF + ".gz")) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.VCF);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.VCF);
                        } else {
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = guessAtUnknowns(myOpenFiles[i]);
                        }
                    } else {
                        alreadyLoaded.add(myOpenFiles[i]);
                        tds = processDatum(myOpenFiles[i], myFileType);
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    StringBuilder builder = new StringBuilder();
                    builder.append("Error loading: ");
                    builder.append(myOpenFiles[i]);
                    builder.append("\n");
                    builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
                    String str = builder.toString();
                    if (isInteractive()) {
                        DialogUtils.showError(str, getParentFrame());
                    } else {
                        myLogger.error(str);
                    }
                }

                if (tds != null) {
                    result.add(tds);
                    fireDataSetReturned(new PluginEvent(tds, FileLoadPlugin.class));
                }

            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public DataSet guessAtUnknowns(String filename) {

        TasselFileType guess = TasselFileType.Sequence;
        DataSet tds = null;

        try {
            BufferedReader br = null;
            if (filename.startsWith("http")) {
                URL url = new URL(filename);
                br = new BufferedReader(new InputStreamReader(url.openStream()));
            } else {
                br = new BufferedReader(new FileReader(filename));
            }

            String line1 = br.readLine().trim();
            String[] sval1 = line1.split("\\s");
            String line2 = br.readLine().trim();
            String[] sval2 = line2.split("\\s");
            boolean lociMatchNumber = false;
            if (!sval1[0].startsWith("<") && (sval1.length == 2) && (line1.indexOf(':') < 0)) {
                int countLoci = Integer.parseInt(sval1[1]);
                if (countLoci == sval2.length) {
                    lociMatchNumber = true;
                }
            }
            if (sval1[0].equalsIgnoreCase("<Annotated>")) {
                guess = TasselFileType.Annotated;
            } else if (line1.startsWith("<") || line1.startsWith("#")) {
                boolean isTrait = false;
                boolean isMarker = false;
                boolean isNumeric = false;
                boolean isMap = false;
                Pattern tagPattern = Pattern.compile("[<>\\s]+");
                String[] info1 = tagPattern.split(line1);
                String[] info2 = tagPattern.split(line2);
                if (info1.length > 1) {
                    if (info1[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info1[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info1[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    } else if (info1[1].toUpperCase().startsWith("MAP")) {
                        isMap = true;
                    }
                }
                if (info2.length > 1) {
                    if (info2[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info2[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info2[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    } else if (info2[1].toUpperCase().startsWith("MAP")) {
                        isMap = true;
                    }
                } else {
                    guess = null;
                    String inline = br.readLine();
                    while (guess == null && inline != null && (inline.startsWith("#") || inline.startsWith("<"))) {
                        if (inline.startsWith("<")) {
                            String[] info = tagPattern.split(inline);
                            if (info[1].toUpperCase().startsWith("MARKER")) {
                                isMarker = true;
                            } else if (info[1].toUpperCase().startsWith("TRAIT")) {
                                isTrait = true;
                            } else if (info[1].toUpperCase().startsWith("NUMER")) {
                                isNumeric = true;
                            } else if (info[1].toUpperCase().startsWith("MAP")) {
                                isMap = true;
                            }
                        }
                    }
                }
                if (isTrait || (isMarker && isNumeric)) {
                    guess = TasselFileType.Phenotype;
                } else if (isMarker) {
                    guess = TasselFileType.Polymorphism;
                } else if (isMap) {
                    guess = TasselFileType.GeneticMap;
                } else {
                    throw new IOException("Improperly formatted header. Data will not be imported.");
                }
            } else if ((line1.startsWith(">")) || (line1.startsWith(";"))) {
                guess = TasselFileType.Fasta;
            } else if (sval1.length == 1) {
                guess = TasselFileType.SqrMatrix;
            } else if (line1.indexOf(':') > 0) {
                guess = TasselFileType.Polymorphism;
            } else if ((sval1.length == 2) && (lociMatchNumber)) {
                guess = TasselFileType.Polymorphism;
            } else if ((line1.startsWith("#Nexus")) || (line1.startsWith("#NEXUS")) || (line1.startsWith("CLUSTAL"))
                    || ((sval1.length == 2) && (sval2.length == 2))) {
                guess = TasselFileType.Sequence;
            } else if (sval1.length == 3) {
                guess = TasselFileType.Numerical;
            }

            myLogger.info("guessAtUnknowns: type: " + guess);
            tds = processDatum(filename, guess);

            br.close();
        } catch (Exception e) {
        }

        return tds;

    }

    private DataSet processDatum(String inFile, TasselFileType theFT) {
        Object result = null;
        String suffix = null;
        try {
            switch (theFT) {
                case Hapmap: {
                    suffix = FILE_EXT_HAPMAP;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_HAPMAP_GZ;
                    }
                    result = ImportUtils.readFromHapmap(inFile, myIsSBit, this);
                    break;
                }
                case HDF5: {
                    suffix = FILE_EXT_HDF5;
//                    result = BitAlignmentHDF5.getInstance(inFile);
                    result = ImportUtils.readGuessFormat(inFile, true);
                    break;
                }
                case VCF: {
                    suffix = FILE_EXT_VCF;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_VCF + ".gz";
                    }
                    result = ImportUtils.readFromVCF(inFile, this);
                    break;
                }
                case Sequence: {
                    result = ReadSequenceAlignmentUtils.readBasicAlignments(inFile, 40);
                    break;
                }
                case Polymorphism: {
                    result = ReadPolymorphismUtils.readPolymorphismFile(inFile);
                    break;
                }
                case Annotated: {
                    result = ReadPhenotypeUtils.readAnnotatedAlignment(inFile);
                    break;
                }
                case Numerical: {
                    result = ReadPhenotypeUtils.readNumericalAlignment(inFile);
                    break;
                }
                case Fasta: {
                    result = ImportUtils.readFasta(inFile, myIsSBit);
                    break;
                }
                case SqrMatrix: {
                    result = ReadDistanceMatrix.readDistanceMatrix(inFile);
                    break;
                }
                case Phenotype: {
                    result = ReadPhenotypeUtils.readGenericFile(inFile);
                    break;
                }
                case GeneticMap: {
                    result = ReadPolymorphismUtils.readGeneticMapFile(inFile);
                    break;
                }
                case Serial: {
                    result = ImportUtils.readAlignmentFromSerialGZ(inFile);
                    break;
                }
                case Table: {
                    result = TableReportUtils.readDelimitedTableReport(inFile, "\t");
                    break;
                }
            }
        } catch (Exception e) {

            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Error loading: ");
            builder.append(inFile);
            builder.append("\n");
            builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
            String str = builder.toString();
            if (isInteractive()) {
                DialogUtils.showError(str, getParentFrame());
            } else {
                myLogger.error(str);
            }

        }
        if (result != null) {
            String theComment = "";
            if (result instanceof Report) {
                StringWriter sw = new StringWriter();
                ((Report) result).report(new PrintWriter(sw));
                theComment = sw.toString();
            }

            String name = Utils.getFilename(inFile, suffix);

            Datum td = new Datum(name, result, theComment);
            //todo need to add logic of directories.
            DataSet tds = new DataSet(td, this);
            return tds;
        }
        return null;
    }

    /**
     * Provides a open filer that remember the last location something was
     * opened from
     */
    private File[] getOpenFilesByChooser() {
        myOpenFileChooser.setMultiSelectionEnabled(true);
        File[] lopenFiles = null;
        myOpenFileChooser.setVisible(true);
        int returnVal = myOpenFileChooser.showOpenDialog(getParentFrame());
        if (returnVal == JFileChooser.OPEN_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            lopenFiles = myOpenFileChooser.getSelectedFiles();
            TasselPrefs.putOpenDir(myOpenFileChooser.getCurrentDirectory().getPath());
        }
        return lopenFiles;
    }

    public String[] getOpenFiles() {
        return myOpenFiles;
    }

    public void setOpenFiles(File[] openFiles) {

        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
            return;
        }

        myOpenFiles = new String[openFiles.length];
        for (int i = 0; i < openFiles.length; i++) {
            myOpenFiles[i] = openFiles[i].getPath();
        }

    }

    public void setOpenFiles(String[] openFiles) {
        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
        } else {
            myOpenFiles = openFiles;
        }
    }

    public TasselFileType getTheFileType() {
        return myFileType;
    }

    public void setTheFileType(TasselFileType theFileType) {
        myFileType = theFileType;
    }

    public boolean getIsFileCreatedSBit() {
        return myIsSBit;
    }

    public void setIsFileCreatedSBit(boolean value) {
        myIsSBit = value;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FileLoadPlugin.class.getResource("images/LoadFile.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Load";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load data from files on your computer.";
    }
}

/**
 * <p>Title: TASSEL</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p>Company: USDA-ARS</p>
 *
 * @author Edward Buckler
 * @version 1.0
 */
class FileLoadPluginDialog extends JDialog {

    boolean isCancel = true;
    ButtonGroup conversionButtonGroup = new ButtonGroup();
    JRadioButton hapMapRadioButton = new JRadioButton("Load Hapmap");
    JRadioButton hdf5RadioButton = new JRadioButton("Load HDF5");
    JRadioButton vcfRadioButton = new JRadioButton("Load VCF");
    JRadioButton plinkRadioButton = new JRadioButton("Load Plink");
    JRadioButton sequenceAlignRadioButton = new JRadioButton("Load sequence alignment (phylip, NEXUS)");
    JRadioButton fastaRadioButton = new JRadioButton("Load FASTA file");
    JRadioButton polymorphismAlignRadioButton = new JRadioButton("Load polymorphism alignment (custom)");
    JRadioButton annotatedAlignRadioButton = new JRadioButton("Load annotated alignment (custom)");
    JRadioButton numericalRadioButton = new JRadioButton("Load numerical trait data or covariates");
    JRadioButton loadMatrixRadioButton = new JRadioButton("Load square numerical matrix (eg. kinship) (phylip)");
    JRadioButton guessRadioButton = new JRadioButton("I will make my best guess and try.");
    JRadioButton flapjackRadioButton = new JRadioButton("Load Flapjack");
    JRadioButton geneticMapRadioButton = new JRadioButton("Load a Genetic Map");
    JRadioButton tableReportRadioButton = new JRadioButton("Load a Table Report");

    public FileLoadPluginDialog() {
        super((Frame) null, "File Loader", true);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        setTitle("File Loader");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);


        Container contentPane = getContentPane();

        BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
        contentPane.setLayout(layout);

        JPanel main = getMain();

        contentPane.add(main);

        pack();

        setResizable(false);

        conversionButtonGroup.add(flapjackRadioButton);
        conversionButtonGroup.add(hapMapRadioButton);
        conversionButtonGroup.add(hdf5RadioButton);
        conversionButtonGroup.add(vcfRadioButton);
        conversionButtonGroup.add(plinkRadioButton);
        conversionButtonGroup.add(sequenceAlignRadioButton);
        conversionButtonGroup.add(fastaRadioButton);
        conversionButtonGroup.add(polymorphismAlignRadioButton);
        conversionButtonGroup.add(loadMatrixRadioButton);
        conversionButtonGroup.add(numericalRadioButton);
        conversionButtonGroup.add(annotatedAlignRadioButton);
        conversionButtonGroup.add(geneticMapRadioButton);
        conversionButtonGroup.add(tableReportRadioButton);
        conversionButtonGroup.add(guessRadioButton);
        guessRadioButton.setSelected(true);

    }

    private JPanel getMain() {

        JPanel inputs = new JPanel();
        BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
        inputs.setLayout(layout);
        inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getLabel());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getOptionPanel());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getButtons());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        return inputs;

    }

    private JPanel getLabel() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);

        JLabel jLabel1 = new JLabel("Choose File Type to Load.");
        jLabel1.setFont(new Font("Dialog", Font.BOLD, 18));
        result.add(jLabel1);

        return result;

    }

    private JPanel getOptionPanel() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        result.setBorder(BorderFactory.createEtchedBorder());

        result.add(hapMapRadioButton);
        result.add(hdf5RadioButton);
        result.add(vcfRadioButton);
        result.add(plinkRadioButton);
        result.add(flapjackRadioButton);
        result.add(sequenceAlignRadioButton);
        result.add(fastaRadioButton);
        result.add(polymorphismAlignRadioButton);
        result.add(annotatedAlignRadioButton);
        result.add(numericalRadioButton);
        result.add(loadMatrixRadioButton);
        result.add(geneticMapRadioButton);
        result.add(tableReportRadioButton);
        result.add(guessRadioButton);

        result.add(Box.createRigidArea(new Dimension(1, 20)));

        return result;

    }

    private JPanel getButtons() {

        JButton okButton = new JButton();
        JButton cancelButton = new JButton();

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });

        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });

        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

        result.add(okButton);

        result.add(cancelButton);

        return result;

    }

    public FileLoadPlugin.TasselFileType getTasselFileType() {
        if (hapMapRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Hapmap;
        }
        if (hdf5RadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.HDF5;
        }
        if (vcfRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.VCF;
        }
        if (plinkRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Plink;
        }
        if (flapjackRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Flapjack;
        }
        if (sequenceAlignRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Sequence;
        }
        if (fastaRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Fasta;
        }
        if (polymorphismAlignRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Polymorphism;
        }
        if (loadMatrixRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.SqrMatrix;
        }
        if (numericalRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Unknown;
        }
        if (annotatedAlignRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Annotated;
        }
        if (geneticMapRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.GeneticMap;
        }
        if (tableReportRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Table;
        }
        return FileLoadPlugin.TasselFileType.Unknown;
    }

    public void okButton_actionPerformed(ActionEvent e) {
        isCancel = false;
        setVisible(false);
    }

    public void cancelButton_actionPerformed(ActionEvent e) {
        isCancel = true;
        setVisible(false);
    }

    public boolean isCancel() {
        return isCancel;
    }
}
