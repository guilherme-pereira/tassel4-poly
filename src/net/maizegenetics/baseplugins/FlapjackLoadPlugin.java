/*
 * FlapjackLoadPlugin.java
 *
 * Created on December 18, 2009
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.*;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import net.maizegenetics.prefs.TasselPrefs;

import net.maizegenetics.util.Utils;

import javax.swing.*;

import java.awt.Container;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;

import java.io.*;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class FlapjackLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FlapjackLoadPlugin.class);
    private String myGenoFile = null;
    private String myMapFile = null;
    private String myChromosome = null;
    private boolean hasHetSeparator;
    private boolean hasNucleotides = true;
    private boolean isPhysicalMap = true;
    private String hetSeparator;
    private String missingCharacter;

    /** Creates a new instance of FlapjackLoadPlugin */
    public FlapjackLoadPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        if (isInteractive()) {
            FlapjackPluginDialog theDialog = new FlapjackPluginDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }
            hasHetSeparator = theDialog.chkHetsSeparated.isSelected();
            hasNucleotides = theDialog.chkMarkersAreNucleotides.isSelected();
            hetSeparator = theDialog.txtHetsep.getText().trim();
            missingCharacter = theDialog.txtMissing.getText().trim();
            isPhysicalMap = theDialog.radioPhysical.isSelected();
            theDialog.dispose();
        }

        if ((myGenoFile == null) || (myGenoFile.length() == 0)) {
            return null;
        }
        if ((myMapFile == null) || (myMapFile.length() == 0)) {
            return null;
        }

        try {
            DataSet result = loadFile(myGenoFile, myMapFile, null);
            return result;
        } catch (Exception e) {
            String msg = "Flapjack files " + myGenoFile + " and " + myMapFile + " failed to load. "
                    + "Make sure the import options are properly set.";
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), msg, "Error uploading Flapjack files", JOptionPane.ERROR_MESSAGE);
            } else {
                myLogger.error(msg);
            }
            return null;
        }

    }

    public String getGenoFile() {
        return myGenoFile;
    }

    public String getMapFile() {
        return myMapFile;
    }

    public String getChromosome() {
        return myChromosome;
    }

    public void setGenoFile(String filename) {
        myGenoFile = filename;
    }

    public void setMapFile(String filename) {
        myMapFile = filename;
    }

    public void setChromosome(String chromosome) {
        myChromosome = chromosome;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Load Flapjack";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load Flapjack Files";
    }

    public DataSet loadFile(String theGenoFile, String theMapFile, String chromosome) {

        Alignment result = null;
        //Terry - fix this
        try {
            if (isPhysicalMap) {
                //result = ImportUtils.readFromFlapjackPhysical(theGenoFile, theMapFile, hasHetSeparator, hasNucleotides, missingCharacter, hetSeparator);
            } else {
                //result = ImportUtils.readFromFlapjackGenetic(theGenoFile, theMapFile, hasHetSeparator, hasNucleotides, missingCharacter, hetSeparator);
            }
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        Datum td = new Datum(Utils.getFilename(theGenoFile, FileLoadPlugin.FILE_EXT_FLAPJACK_GENO), result, null);
        DataSet tds = new DataSet(td, this);
        fireDataSetReturned(new PluginEvent(tds, FlapjackLoadPlugin.class));

        return tds;

    }

    public void hasHetSeparator(boolean hasSeparator) {
        hasHetSeparator = hasSeparator;
    }

    public void hasNucleotides(boolean nucleotides) {
        hasNucleotides = nucleotides;
    }

    public void setHetSeparator(String separator) {
        hetSeparator = separator;
    }

    public void setMissingCharacter(String missing) {
        missingCharacter = missing;
    }

    public void hasPhysicalMap(boolean physicalmap) {
        isPhysicalMap = physicalmap;
    }

    class FlapjackPluginDialog extends JDialog {

        private JPanel main = null;
        private final static int TEXT_FIELD_WIDTH = 30;
        private JTextField myMapFileField = null;
        private JTextField myGenoFileField = null;
        private JTextField myChromosomeField = null;
        private boolean myIsCancel = false;
        private JButton myMapFileBrowseButton = null;
        private JButton myGenoFileBrowseButton = null;
        private final JCheckBox chkHetsSeparated = new JCheckBox("Heterozygotes are separated by a character (eg, A/T rather than AT)", true);
        private final JCheckBox chkMarkersAreNucleotides = new JCheckBox("Markers are nucleotides (ie, A, C, G, or T)", true);
        private final JTextField txtHetsep = new JTextField(5);
        private final JTextField txtMissing = new JTextField(5);
        private final JRadioButton radioPhysical = new JRadioButton("Physical", true);
        private final JRadioButton radioGenetic = new JRadioButton("Genetic", false);

        public FlapjackPluginDialog() {
            super((Frame) null, "File Loader", true);
            try {
                createDialog();
                pack();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        private void createDialog() {

            myGenoFileField = new JTextField(TEXT_FIELD_WIDTH);
            myGenoFileField.setText("(Select .geno File)");

            myMapFileField = new JTextField(TEXT_FIELD_WIDTH);
            myMapFileField.setText("(Select .MAP File)");

            myGenoFileBrowseButton = new JButton("Browse...");
            myMapFileBrowseButton = new JButton("Browse...");
            txtHetsep.setText("/");
            txtMissing.setText("-");

            Frame frame = JOptionPane.getFrameForComponent(getParentFrame());
            setLocationRelativeTo(frame);

            setTitle("Load Flapjack");
            setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
            setUndecorated(false);
            getRootPane().setWindowDecorationStyle(JRootPane.NONE);

            Container contentPane = getContentPane();

            contentPane.setLayout(new GridBagLayout());

            addComponentListener(new ComponentAdapter() {

                public void componentShown(ComponentEvent ce) {
                    myGenoFileField.requestFocusInWindow();
                }
            });

            myGenoFileBrowseButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    JFileChooser myFileChooser = new JFileChooser(TasselPrefs.getOpenDir());
                    myFileChooser.setDialogTitle("Open a Flapjack geno file");
                    if (myFileChooser.showOpenDialog(main) == JFileChooser.APPROVE_OPTION) {
                        File file = myFileChooser.getSelectedFile();
                        myGenoFileField.setText(file.getPath());
                        TasselPrefs.putOpenDir(myFileChooser.getCurrentDirectory().getPath());
                    }
                }
            });

            myMapFileBrowseButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    JFileChooser myFileChooser = new JFileChooser(TasselPrefs.getOpenDir());
                    myFileChooser.setDialogTitle("Open a Flapjack map file");
                    if (myFileChooser.showOpenDialog(main) == JFileChooser.APPROVE_OPTION) {
                        File file = myFileChooser.getSelectedFile();
                        myMapFileField.setText(file.getPath());
                        TasselPrefs.putOpenDir(myFileChooser.getCurrentDirectory().getPath());
                    }
                }
            });

            GridBagConstraints gbc = new GridBagConstraints();

            gbc.gridx = 0;
            gbc.gridy = 0;
            gbc.insets = new Insets(6, 15, 4, 4);
            gbc.anchor = GridBagConstraints.WEST;
            contentPane.add(new JLabel("Genotype File"), gbc);
            gbc.gridx++;
            gbc.gridwidth = 2;
            gbc.insets.left = 4;
            gbc.anchor = GridBagConstraints.CENTER;
            contentPane.add(myGenoFileField, gbc);
            gbc.gridx += 2;
            gbc.gridwidth = 1;
            gbc.insets.right = 10;
            gbc.anchor = GridBagConstraints.WEST;
            contentPane.add(myGenoFileBrowseButton, gbc);

            gbc.gridy++;
            gbc.gridx = 0;
            gbc.anchor = GridBagConstraints.WEST;
            gbc.insets.left = 15;
            contentPane.add(new JLabel("Map File"), gbc);
            gbc.gridx++;
            gbc.gridwidth = 2;
            gbc.insets.left = 4;
            gbc.anchor = GridBagConstraints.CENTER;
            contentPane.add(myMapFileField, gbc);
            gbc.gridx += 2;
            gbc.gridwidth = 1;
            gbc.insets.right = 10;
            gbc.anchor = GridBagConstraints.WEST;
            contentPane.add(myMapFileBrowseButton, gbc);

            gbc.gridy++;
            gbc.gridx = 0;
            gbc.gridwidth = 3;
            gbc.anchor = GridBagConstraints.WEST;
            gbc.insets.left = 15;
            gbc.insets.top = 25;
            contentPane.add(chkHetsSeparated, gbc);
            gbc.gridy++;
            gbc.insets.top = 5;
            contentPane.add(chkMarkersAreNucleotides, gbc);

            gbc.gridwidth = 2;
            gbc.gridy++;
            contentPane.add(new JLabel("Heterozygote separator character"), gbc);
            gbc.gridx += 2;
            gbc.gridwidth = 1;
            gbc.insets.left = 5;
            contentPane.add(txtHetsep, gbc);
            gbc.gridy++;
            gbc.gridx = 0;
            gbc.gridwidth = 2;
            gbc.insets.left = 15;
            contentPane.add(new JLabel("Missing data character"), gbc);
            gbc.gridx += 2;
            gbc.insets.left = 5;
            gbc.gridwidth = 1;
            contentPane.add(txtMissing, gbc);

            //radio buttons for map type
            ButtonGroup bg = new ButtonGroup();
            bg.add(radioPhysical);
            bg.add(radioGenetic);
            JPanel mapTypePanel = new JPanel();
            mapTypePanel.setBorder(BorderFactory.createTitledBorder("Map Type"));
            mapTypePanel.setLayout(new BoxLayout(mapTypePanel, BoxLayout.Y_AXIS));
            mapTypePanel.add(radioPhysical);
            mapTypePanel.add(radioGenetic);
            gbc.gridy++;
            gbc.gridx = 0;
            gbc.gridwidth = 5;
            gbc.insets.top = 25;
            gbc.fill = GridBagConstraints.HORIZONTAL;
            gbc.anchor = GridBagConstraints.WEST;
            contentPane.add(mapTypePanel, gbc);

            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.X_AXIS));
            buttonPanel.add(getOkButton());
            buttonPanel.add(Box.createHorizontalStrut(30));
            buttonPanel.add(getCancelButton());

            gbc.gridy++;
            gbc.gridx = 0;
            gbc.gridwidth = 6;
            gbc.insets.top = 25;
            gbc.anchor = GridBagConstraints.CENTER;
            gbc.fill = GridBagConstraints.NONE;
            contentPane.add(buttonPanel, gbc);

            pack();
        }

        private JButton getResetButton() {

            JButton resetButton = new JButton("Reset");
            resetButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    myGenoFileField.setText("");
                    myMapFileField.setText("");
                    myChromosomeField.setText("");
                }
            });

            return resetButton;

        }

        private JButton getCancelButton() {

            JButton cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    myIsCancel = true;
                    setVisible(false);
                }
            });

            return cancelButton;

        }

        private JButton getOkButton() {

            JButton okButton = new JButton();
            okButton.setText("Import");
            okButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    myGenoFile = myGenoFileField.getText();
                    myMapFile = myMapFileField.getText();
//                    myChromosome = myChromosomeField.getText();
                    myIsCancel = false;
                    setVisible(false);
                }
            });

            return okButton;

        }

        public boolean isCancel() {
            return myIsCancel;
        }
    }
}
