/*
 * PlinkLoadPlugin.java
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
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import java.io.*;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class PlinkLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(PlinkLoadPlugin.class);
    private String myPedFile = null;
    private String myMapFile = null;
    private String myChromosome = null;

    /** Creates a new instance of PlinkLoadPlugin */
    public PlinkLoadPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        if (isInteractive()) {
            PlinkLoadPluginDialog theDialog = new PlinkLoadPluginDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }

            theDialog.dispose();
        }

        if ((myPedFile == null) || (myPedFile.length() == 0)) {
            return null;
        }
        if ((myMapFile == null) || (myMapFile.length() == 0)) {
            return null;
        }

        DataSet result = loadFile(myPedFile, myMapFile, myChromosome);

        return result;

    }

    public String getPedFile() {
        return myPedFile;
    }

    public String getMapFile() {
        return myMapFile;
    }

    public String getChromosome() {
        return myChromosome;
    }

    public void setPedFile(String filename) {
        myPedFile = filename;
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
        return "Load Plink";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load Plink Files";
    }

    public DataSet loadFile(String thePedFile, String theMapFile, String chromosome) {

        // Terry - fix this.
        Alignment result = null;
        //Alignment result = ImportUtils.readFromPLINK(thePedFile, theMapFile);
        Datum td = new Datum(Utils.getFilename(thePedFile, FileLoadPlugin.FILE_EXT_PLINK_PED), result, null);
        DataSet tds = new DataSet(td, this);
        fireDataSetReturned(new PluginEvent(tds, PlinkLoadPlugin.class));

        return tds;

    }

    class PlinkLoadPluginDialog extends JDialog {

        private JPanel main = null;
        private final static int TEXT_FIELD_WIDTH = 30;
        private JTextField myMapFileField = null;
        private JTextField myPedFileField = null;
        private JTextField myChromosomeField = null;
        private boolean myIsCancel = false;

        private JButton myMapFileBrowseButton = null;
        private JButton myPedFileBrowseButton = null;

        private JFileChooser myFileChooser = null;

        public PlinkLoadPluginDialog() {
            super((Frame) null, "File Loader", true);
            try {
                createDialog();
                pack();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        private void createDialog() {

            myPedFileField = new JTextField(TEXT_FIELD_WIDTH);
            myPedFileField.setText("(Select .PED File)");

            myMapFileField = new JTextField(TEXT_FIELD_WIDTH);
            myMapFileField.setText("(Select .MAP File)");
            
//            myChromosomeField = new JTextField(TEXT_FIELD_WIDTH);
            myFileChooser = new JFileChooser(TasselPrefs.getOpenDir());

            myPedFileBrowseButton = new JButton("Browse...");
            myMapFileBrowseButton = new JButton("Browse...");

            myPedFileField.addFocusListener(new FocusAdapter() {

                public void focusLost(FocusEvent e) {
                    String temp = myMapFileField.getText();
                    if ((temp == null) || (temp.length() == 0)) {
                        myMapFileField.setText(myPedFileField.getText());
                    }
                }
            });

            Frame frame = JOptionPane.getFrameForComponent(getParentFrame());
            setLocationRelativeTo(frame);

            setTitle("Load Plink");
            setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
            setUndecorated(false);
            getRootPane().setWindowDecorationStyle(JRootPane.NONE);

            Container contentPane = getContentPane();

            BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
            contentPane.setLayout(layout);

            main = getMain();

            contentPane.add(main);

            pack();

            setResizable(false);

            addComponentListener(new ComponentAdapter() {

                public void componentShown(ComponentEvent ce) {
                    myPedFileField.requestFocusInWindow();
                }
            });

            myPedFileBrowseButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    if (myFileChooser.showOpenDialog(main) == JFileChooser.APPROVE_OPTION) {
                        File file = myFileChooser.getSelectedFile();
                        myPedFileField.setText(file.getPath());
                        TasselPrefs.putOpenDir(myFileChooser.getCurrentDirectory().getPath());
                    }
                }

            });

            myMapFileBrowseButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    if (myFileChooser.showOpenDialog(main) == JFileChooser.APPROVE_OPTION) {
                        File file = myFileChooser.getSelectedFile();
                        myMapFileField.setText(file.getPath());
                        TasselPrefs.putOpenDir(myFileChooser.getCurrentDirectory().getPath());
                    }
                }

            });
        }

        private JPanel getMain() {

            JPanel inputs = new JPanel();
            BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
            inputs.setLayout(layout);
            inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);

            inputs.add(getLine("Ped File:", myPedFileField, myPedFileBrowseButton));
            inputs.add(getLine("Map File:", myMapFileField, myMapFileBrowseButton));
//            inputs.add(getLine("Chromosome:", myChromosomeField));

            inputs.add(Box.createRigidArea(new Dimension(1, 20)));

            inputs.add(getButtons());

            return inputs;

        }

        private JButton getResetButton() {

            JButton resetButton = new JButton("Reset");
            resetButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    myPedFileField.setText("");
                    myMapFileField.setText("");
                    myChromosomeField.setText("");
                }
            });

            return resetButton;

        }

        private JPanel getLine(String label, JTextField ref) {

            JPanel result = new JPanel(new FlowLayout(FlowLayout.LEFT));

            result.add(new JLabel(label));
            ref.setEditable(true);
            ref.setHorizontalAlignment(JTextField.LEFT);
            ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
            ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
            ref.setMaximumSize(ref.getPreferredSize());
            result.add(ref);

            return result;

        }

        private JPanel getLine(String label, JTextField ref, JButton button) {

            JPanel result = new JPanel(new FlowLayout(FlowLayout.LEFT));

            result.add(new JLabel(label));
            ref.setEditable(true);
            ref.setHorizontalAlignment(JTextField.LEFT);
            ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
            ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
            ref.setMaximumSize(ref.getPreferredSize());
            result.add(ref);
            result.add(button);

            return result;

        }

        private JPanel getButtons() {

            JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

            result.add(getOkButton());

            result.add(getCancelButton());

            result.add(getResetButton());

            return result;

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
            okButton.setText("Ok");
            okButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    myPedFile = myPedFileField.getText();
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
