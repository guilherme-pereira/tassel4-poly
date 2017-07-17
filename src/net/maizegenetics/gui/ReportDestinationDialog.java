package net.maizegenetics.gui;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

public class ReportDestinationDialog extends JDialog {

    private JCheckBox chkFileOutput; //checked if the output is to be sent to a file
    private JCheckBox chkRestrictOutput; //checked if the ouput is to be filtered on maxP
    private JTextField txtFilename;	//a text field holding the output file name
    private JTextField txtMaxP; //a text field holding the max P value to be output
    private double maxp = 1;  //maximum pvalue to be returned
    private JDialog thisDialog = this;
    private boolean okayChecked = false;

    //for testing
    public static void main(String[] args) {
        ReportDestinationDialog rdd = new ReportDestinationDialog();
        rdd.setVisible(true);
        System.exit(-1);
    }

    /**
     * Creates a dialog that allows the user to send results to a file rather than to save them to memory.
     * Very large output sets can fill up memory and cause TASSEL to crash otherwise.
     *
     */
    public ReportDestinationDialog() {
        this.setModal(true);
        this.setTitle("Choose Output Format");
        String msg = "By default all output will be stored in memory. "
                + "If the output is large enough it could exceed the memory space available to TASSEL and cause TASSEL to crash. "
                + "If you think this could be a problem, you can choose to write the output to a text file instead. "
                + "Enter only the file base name. The dataset name and file type will be appended when the report is saved.";
        Container contentPane = this.getContentPane();
        contentPane.setLayout(new GridBagLayout());

        txtFilename = new JTextField(30);
        txtMaxP = new JTextField(8);
        txtMaxP.setText("1e-3");
        final JTextArea txtWarning = new JTextArea(msg);
        txtWarning.setLineWrap(true);
        txtWarning.setWrapStyleWord(true);
        txtWarning.setPreferredSize(new Dimension(500, 150));
        txtWarning.setEditable(false);
        final JLabel lblFilename = new JLabel("Output base file name:");
        final JLabel lblMaxp = new JLabel("Do not keep p-values larger than");
        final JButton btnBrowse = new JButton("Browse");
        final JButton btnOkay = new JButton("Okay");

        lblFilename.setEnabled(false);
        lblMaxp.setEnabled(false);
        btnBrowse.setEnabled(false);
        txtFilename.setEnabled(false);
        txtMaxP.setEnabled(false);

        chkFileOutput = new JCheckBox("Write output to file", false);
        chkFileOutput.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent arg0) {
                if (chkFileOutput.isSelected()) {
                    lblFilename.setEnabled(true);
                    txtFilename.setEnabled(true);
                    btnBrowse.setEnabled(true);
                } else {
                    lblFilename.setEnabled(false);
                    txtFilename.setEnabled(false);
                    btnBrowse.setEnabled(false);
                }
            }
        });

        chkRestrictOutput = new JCheckBox("Filter output on p-value", false);
        chkRestrictOutput.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                if (chkRestrictOutput.isSelected()) {
                    lblMaxp.setEnabled(true);
                    txtMaxP.setEnabled(true);
                } else {
                    lblMaxp.setEnabled(false);
                    txtMaxP.setEnabled(false);
                }
            }
        });

        btnBrowse.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent evt) {
                JFileChooser chooser = new JFileChooser();
                chooser.setDialogTitle("Choose a File for Output");
                int action = chooser.showSaveDialog(thisDialog);

                if (action == JFileChooser.APPROVE_OPTION) {
                    File theFile = chooser.getSelectedFile();
                    txtFilename.setText(theFile.getPath());
                }
            }
        });

        btnOkay.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent evt) {
                //if pmax is checked and the value in the box is not a valid number generate an error
                boolean validp = true;
                okayChecked = true;
                if (chkRestrictOutput.isSelected()) {
                    try {
                        maxp = Double.parseDouble(txtMaxP.getText());
                    } catch (Exception e) {
                        validp = false;
                    }
                }

                if (validp) {
                    setVisible(false);
                } else {
                    JOptionPane.showMessageDialog(thisDialog, "Max P is not a valid number.");
                }
            }
        });

        GridBagConstraints gbc = new GridBagConstraints();
        int y = 0;
        gbc.gridx = 0;
        gbc.gridy = y++;
        gbc.gridwidth = 1;
        gbc.insets = new Insets(20, 5, 20, 5);
        contentPane.add(txtWarning, gbc);

        JPanel pnlFile = new JPanel(new GridBagLayout());
        pnlFile.setBorder(BorderFactory.createLineBorder(Color.BLACK));

        gbc.gridx = 0;
        gbc.gridy = y++;
        gbc.insets = new Insets(10, 10, 10, 10);
        contentPane.add(pnlFile, gbc);

        gbc.gridy = 0;
        gbc.gridx = 0;
        gbc.gridwidth = 3;
        gbc.insets = new Insets(15, 5, 5, 5);
        pnlFile.add(chkFileOutput, gbc);

        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 1;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlFile.add(lblFilename, gbc);

        gbc.gridx = 1;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlFile.add(txtFilename, gbc);

        gbc.gridx = 2;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlFile.add(btnBrowse, gbc);

        JPanel pnlMaxp = new JPanel(new GridBagLayout());
        pnlMaxp.setBorder(BorderFactory.createLineBorder(Color.BLACK));

        gbc.gridx = 0;
        gbc.gridy = y++;
        gbc.gridwidth = 1;
        gbc.insets = new Insets(10, 10, 10, 10);
        contentPane.add(pnlMaxp, gbc);

        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 3;
        gbc.insets = new Insets(15, 5, 5, 5);
        pnlMaxp.add(chkRestrictOutput, gbc);

        gbc.gridx = 1;
        gbc.gridy = y++;
        gbc.gridwidth = 1;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlMaxp.add(lblMaxp, gbc);

        gbc.gridx = 2;
        gbc.insets = new Insets(5, 5, 5, 5);
        pnlMaxp.add(txtMaxP, gbc);

        gbc.gridx = 0;
        gbc.gridy = y++;
        gbc.gridwidth = 1;
        gbc.insets = new Insets(15, 5, 20, 5);
        contentPane.add(btnOkay, gbc);

        this.pack();
    }

    public double getMaxP() {
        return maxp;
    }

    public String getOutputFileName() {
        return txtFilename.getText().trim();
    }

    public boolean wasUseFileChecked() {
        return chkFileOutput.isSelected();
    }

    public boolean wasRestrictOutputChecked() {
        return chkRestrictOutput.isSelected();
    }

    public boolean isOkayChecked() {
        return okayChecked;
    }
}
