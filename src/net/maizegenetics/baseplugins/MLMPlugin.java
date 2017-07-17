package net.maizegenetics.baseplugins;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.pal.alignment.GeneticMap;
import net.maizegenetics.pal.alignment.MarkerPhenotype;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.stats.MLM.*;
import net.maizegenetics.gui.ReportDestinationDialog;

import javax.swing.*;
import java.net.URL;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;



import org.apache.log4j.Logger;

/**
 * Author: Peter Bradbury
 * Date: Jan 17, 2010
 */
public class MLMPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(MLMPlugin.class);
    protected DistanceMatrix kinshipMatrix;
    protected boolean analyzeByColumn;
    protected boolean useP3D = true;
    protected CompressionType compressionType = CompressionType.Optimum;
    protected double compression = 1;
    private boolean writeOutputToFile = false;
    private String outputName = null;
    private boolean filterOutput = false;
    private double maxp = 1;

    public enum CompressionType {

        Optimum, Custom, None
    };

    public MLMPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);

    }

    public DataSet performFunction(DataSet input) {

        try {
            //import a genetic map, if there is one
            GeneticMap myMap;
            List<Datum> maps = input.getDataOfType(GeneticMap.class);
            if (maps.size() > 0) {
                myMap = (GeneticMap) maps.get(0).getData();
            } else {
                myMap = null;
            }

            java.util.List<Datum> alignInList = input.getDataOfType(MarkerPhenotype.class);
            if (alignInList.size() == 0) {
                alignInList = input.getDataOfType(Phenotype.class);
            }
            java.util.List<Datum> kinshipList = input.getDataOfType(DistanceMatrix.class);

            if (alignInList.size() < 1) {
                String message = "Invalid selection. Please select sequence alignment or marker and trait data.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error("performFunction: " + message);
                }
                return null;
            }
            if (kinshipList.size() != 1) {
                String message = "Please select exactly one kinship matrix.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error("performFunction: " + message);
                }
                return null;
            }

            kinshipMatrix = (DistanceMatrix) kinshipList.get(0).getData();
            Iterator<Datum> itr = alignInList.iterator();

            if (isInteractive()) {
                MLMOptionDialog theOptions = new MLMOptionDialog(getParentFrame());

                if (theOptions.runClicked) {
                    useP3D = theOptions.useP3D();
                    compressionType = theOptions.getCompressionType();
                    compression = theOptions.getCompressionLevel();
                    theOptions.dispose();

                    // give the user the option of sending the output to a file
                    ReportDestinationDialog rdd = new ReportDestinationDialog();
                    rdd.setLocationRelativeTo(getParentFrame());
                    rdd.setVisible(true);
                    if (!rdd.isOkayChecked()) {
                        return null;
                    }
                    writeOutputToFile = rdd.wasUseFileChecked();
                    if (writeOutputToFile) {
                        outputName = rdd.getOutputFileName();
                    }
                    filterOutput = rdd.wasRestrictOutputChecked();
                    if (filterOutput) {
                        maxp = rdd.getMaxP();
                    }

                } else {
                    theOptions.dispose();
                    return null;
                }

            } else {         //non-interactive stuff
            }

            List<DataSet> result = new ArrayList<DataSet>();
            while (itr.hasNext()) {

                Datum current = itr.next();

                DataSet tds = null;

                try {
                    if (useP3D) {
                        if (compressionType.equals(CompressionType.Optimum)) {
                            CompressedMLMusingDoubleMatrix theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, true, Double.NaN, myMap);
                            tds = new DataSet(theAnalysis.solve(), this);
                        } else if (compressionType.equals(CompressionType.Custom)) {
                            CompressedMLMusingDoubleMatrix theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, true, compression, myMap);
                            tds = new DataSet(theAnalysis.solve(), this);
                        } else {
                            CompressedMLMusingDoubleMatrix theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, false, true, Double.NaN, myMap);
                            tds = new DataSet(theAnalysis.solve(), this);
                        }
                    } else {
                        if (compressionType.equals(CompressionType.Optimum)) {
                            CompressedMLMusingDoubleMatrix theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, false, Double.NaN, myMap);
                            tds = new DataSet(theAnalysis.solve(), this);
                        } else if (compressionType.equals(CompressionType.Custom)) {
                            CompressedMLMusingDoubleMatrix theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, true, false, compression, myMap);
                            tds = new DataSet(theAnalysis.solve(), this);
                        } else {
                            CompressedMLMusingDoubleMatrix theAnalysis = new CompressedMLMusingDoubleMatrix(this, current, kinshipMatrix, false, false, Double.NaN, myMap);
                            tds = new DataSet(theAnalysis.solve(), this);
                        }
                    }
                } catch (Exception e) {
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
                        e.printStackTrace();
                    } else {
                        System.out.println(e.getMessage());
                        e.printStackTrace();
                    }
                }

                if (tds != null) {
                    result.add(tds);
                    fireDataSetReturned(new PluginEvent(tds, MLMPlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public ImageIcon getIcon() {
        URL imageURL = MLMPlugin.class.getResource("images/Mix.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "MLM";
    }

    public String getToolTipText() {
        return "Association analysis using mixed model";
    }

    //a few stub functions to avoid producing errors in existing pipeline code
    public void setAnalyzeByColumn(boolean analyzeByColumn) {
        this.analyzeByColumn = analyzeByColumn;
    }

    public void setMaximumNumOfIteration(int max) {/*does nothing*/

    }

    public void setFinalIterMarker(boolean myFinalIterMarker) {/*does nothing*/

    }

    public void addFactors(int[] factors) {/*does nothing*/

    }

    public void setColumnTypes(String[] types) {/*does nothing*/

    }

    public void addFactors(String[] names) {/*does nothing*/

    }

    public void updateProgress(int progress) {
        if (progress < 0) {
            progress = 0;
        } else if (progress > 100) {
            progress = 100;
        }
        fireProgress(progress);
    }

    public void setVarCompEst(String value) {
        if (value.equalsIgnoreCase("P3D")) {
            useP3D = true;
        } else if (value.equalsIgnoreCase("EachMarker")) {
            useP3D = false;
        } else {
            throw new IllegalArgumentException("MLMPlugin: setVarCompEst: don't know how to handle value: " + value);
        }
    }

    public void setCompressionType(CompressionType type) {
        compressionType = type;
    }

    public boolean isWriteOutputToFile() {
        return writeOutputToFile;
    }

    public void setWriteOutputToFile(boolean writeOutputToFile) {
        this.writeOutputToFile = writeOutputToFile;
    }

    public String getOutputName() {
        return outputName;
    }

    public void setOutputName(String outputName) {
        this.outputName = outputName;
        this.writeOutputToFile = true;
    }

    public boolean isFilterOutput() {
        return filterOutput;
    }

    public void setFilterOutput(boolean filterOutput) {
        this.filterOutput = filterOutput;
    }

    public double getMaxp() {
        return maxp;
    }

    public void setMaxp(double maxp) {
        this.maxp = maxp;
        this.filterOutput = true;
    }

    public double getCustomCompression() {
        return compression;
    }

    public void setCustomCompression(double value) {
        compression = value;
    }
}

class MLMOptionDialog extends JDialog implements ActionListener {

    JRadioButton btnOptimum, btnCustom, btnNoCompression, btnEachMarker, btnP3D;
    ButtonGroup bgCompress, bgVariance;
    JTextField txtCustom;
    boolean runClicked = false;
    boolean useP3D = true;
    MLMPlugin.CompressionType compressionType = MLMPlugin.CompressionType.Optimum;
    JPanel distancePanel;
    
    MLMOptionDialog(Frame parentFrame) {
        super(parentFrame, true);
        final Frame pframe = parentFrame;
        setTitle("MLM Options");
        setSize(new Dimension(350, 300));
        setLocationRelativeTo(pframe);
        Container theContentPane = getContentPane();
        theContentPane.setLayout(new BorderLayout());
        JPanel compressionPanel = new JPanel(new GridBagLayout());
        compressionPanel.setBorder(BorderFactory.createTitledBorder("Compression Level"));

        //the method radio buttons
        btnOptimum = new JRadioButton("Optimum Level", true);
        btnOptimum.setActionCommand("Optimum");
        btnOptimum.addActionListener(this);
        btnCustom = new JRadioButton("Custom Level:", false);
        btnCustom.setActionCommand("Custom");
        btnCustom.addActionListener(this);
        btnNoCompression = new JRadioButton("No Compression", false);
        btnNoCompression.setActionCommand("None");
        btnNoCompression.addActionListener(this);
        btnEachMarker = new JRadioButton("Re-estimate after each marker", false);
        btnEachMarker.setActionCommand("Eachmarker");
        btnEachMarker.addActionListener(this);
        btnP3D = new JRadioButton("P3D (estimate once)", true);
        btnP3D.setActionCommand("P3D");
        btnP3D.addActionListener(this);

        bgCompress = new ButtonGroup();
        bgCompress.add(btnOptimum);
        bgCompress.add(btnCustom);
        bgCompress.add(btnNoCompression);
        bgVariance = new ButtonGroup();
        bgVariance.add(btnEachMarker);
        bgVariance.add(btnP3D);

        txtCustom = new JTextField(5);
        Insets inset1 = new Insets(5, 15, 5, 5);
        Insets inset2 = new Insets(5, 5, 5, 5);
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridwidth = 2;
        gbc.gridy = 0;
        gbc.weightx = 0;
        gbc.weighty = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = inset1; //top, left, bottom, right
        compressionPanel.add(btnOptimum, gbc);

        gbc.gridy++;
        gbc.gridwidth = 1;
        compressionPanel.add(btnCustom, gbc);
        gbc.gridx++;
        gbc.insets = inset2;
        compressionPanel.add(txtCustom, gbc);

        gbc.gridy++;
        gbc.gridx = 0;
        gbc.gridwidth = 2;
        gbc.insets = inset1;
        compressionPanel.add(btnNoCompression, gbc);

        theContentPane.add(compressionPanel, BorderLayout.NORTH);

        JPanel variancePanel = new JPanel(new GridBagLayout());
        variancePanel.setBorder(BorderFactory.createTitledBorder("Variance Component Estimation"));
        gbc.gridy = 0;
        variancePanel.add(btnP3D, gbc);
        gbc.gridy++;
        variancePanel.add(btnEachMarker, gbc);
        theContentPane.add(variancePanel, BorderLayout.CENTER);

        //the help me button
        JButton btnHelpme = new JButton("Help Me Choose");
        final String msg = "For faster analysis, impute marker values before running MLM and use P3D.\n"
                + "With imputed marker values (no missing data), P3D will be very fast but compression will actually increase execution time somewhat.\n"
                + "However, because compression will improve the overall model fit, it should still be used.\n"
                + "If there is missing marker data, compression with P3D will probably be faster than P3D alone.\n"
                + "With small to moderate sized data sets any of the methods should give reasonable performance. \n"
                + "With large data sets consider imputing marker data then using EMMA with both compression and P3D";
        btnHelpme.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent arg0) {
                JOptionPane.showMessageDialog(pframe, msg, "EMMA methods", JOptionPane.INFORMATION_MESSAGE);
            }
        });


        //the run and cancel buttons
        JButton btnRun = new JButton("Run");
        JButton btnCancel = new JButton("Cancel");
        btnRun.setActionCommand("run");
        btnRun.addActionListener(this);
        btnCancel.setActionCommand("cancel");
        btnCancel.addActionListener(this);

        Box buttonBox = Box.createHorizontalBox();
        buttonBox.add(Box.createGlue());
        buttonBox.add(btnRun);
        buttonBox.add(Box.createHorizontalStrut(50));
        buttonBox.add(btnCancel);
        buttonBox.add(Box.createHorizontalStrut(50));
        buttonBox.add(btnHelpme);
        buttonBox.add(Box.createGlue());
        theContentPane.add(buttonBox, BorderLayout.SOUTH);
        this.pack();
        setLocationRelativeTo(getParent());
        this.setVisible(true);

    }

    public boolean useP3D() {
        return useP3D;
    }

    public MLMPlugin.CompressionType getCompressionType() {
        return compressionType;
    }

    double getCompressionLevel() {
        double comp;
        try {
            comp = Double.parseDouble(txtCustom.getText());
        } catch (Exception e) {
            comp = Double.NaN;
        }
        return comp;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equals("run")) {
            runClicked = true;
            this.setVisible(false);
        } else if (e.getActionCommand().equals("cancel")) {
            this.setVisible(false);
        } else if (e.getActionCommand().equals("Optimum")) {
            compressionType = MLMPlugin.CompressionType.Optimum;
        } else if (e.getActionCommand().equals("Custom")) {
            compressionType = MLMPlugin.CompressionType.Custom;
        } else if (e.getActionCommand().equals("None")) {
            compressionType = MLMPlugin.CompressionType.None;
        } else if (e.getActionCommand().equals("Eachmarker")) {
            useP3D = false;
        } else if (e.getActionCommand().equals("P3D")) {
            useP3D = true;
        }

    }
}

class MLMNewOptionDialog extends JDialog implements ActionListener {
	JCheckBox chkP3D = new JCheckBox("P3D (Compute variance estimates once)", true);
	JCheckBox chkCompression = new JCheckBox("Compression", true);
	JCheckBox chkNoMarkers = new JCheckBox("Do not test markers", false);
	JCheckBox chkUPGMA = new JCheckBox("UPGMA", true);
	JCheckBox chkNJ = new JCheckBox("Neighbor Joining", false);
	JCheckBox chkAvg = new JCheckBox("Average", true);
	JCheckBox chkMin = new JCheckBox("Minimum", false);
	JCheckBox chkMax = new JCheckBox("Maximum", false);
	JCheckBox chkMedian = new JCheckBox("Median", false);

	JTextField txtGroupFrom = new JTextField(5);
	JTextField txtGroupTo = new JTextField(5);
	JTextField txtGroupBy = new JTextField(5);
	JTextField txtGroupNumberList = new JTextField(30);
	
	JTextField txtCompFrom = new JTextField(5);
	JTextField txtCompTo = new JTextField(5);
	JTextField txtCompBy = new JTextField(5);
	JTextField txtCompNumberList = new JTextField(30);
	
	JTextField txtCompressionFrom = new JTextField(5);
	JTextField txtCompressionTo = new JTextField(5);
	JTextField txtCompressionLevels = new JTextField(5);
	
	JRadioButton radioGroupRange = new JRadioButton("Range", false);
	JRadioButton radioGroupList = new JRadioButton("List (comma-separated numbers)", false);
	JRadioButton radioCompRange = new JRadioButton("Range", true);
	JRadioButton radioCompList = new JRadioButton("List (comma-separated numbers)", false);
	
    MLMNewOptionDialog(Frame parentFrame) {
        super(parentFrame, true);
        setTitle("MLM Options");
        setSize(new Dimension(350, 300));
        setLocationRelativeTo(getParent());
        Container theContentPane = getContentPane();
        theContentPane.setLayout(new BorderLayout());
        JPanel optionPanel = new JPanel(new GridBagLayout());
        optionPanel.setBorder(BorderFactory.createTitledBorder("MLM Options"));
        JPanel compressionPanel = new JPanel(new GridBagLayout());
        JPanel groupNumberPanel = new JPanel(new GridBagLayout());
        JPanel centerPanel = new JPanel(new BorderLayout());
        JPanel buttonPanel = new JPanel(new GridBagLayout());

        String gnPanelTitle = "Specify Compression Levels or Group Sizes to Test";
        compressionPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(5, 2, 4, 2), BorderFactory.createTitledBorder("Compression Options")));
        groupNumberPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(4, 2, 8, 2), BorderFactory.createTitledBorder(gnPanelTitle)));

        chkCompression.addActionListener(this);
        chkCompression.setActionCommand("compress");
        chkNoMarkers.addActionListener(this);
        chkNoMarkers.setActionCommand("nomarkers");
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.insets = new Insets(2,2,2,2);
        optionPanel.add(chkP3D, gbc);
        gbc.gridy = 1;
        optionPanel.add(chkCompression, gbc);
        gbc.gridy = 2;
        optionPanel.add(chkNoMarkers, gbc);
		
        theContentPane.add(optionPanel, BorderLayout.NORTH);
        
        ButtonGroup numberGroup = new ButtonGroup();
        numberGroup.add(radioGroupRange);
        numberGroup.add(radioGroupList);
        numberGroup.add(radioCompRange);
        numberGroup.add(radioCompList);
        
        gbc.gridy = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(2,4,2,4);
        compressionPanel.add(new JLabel("Clustering Algorithm"), gbc);
        gbc.gridx++;
        compressionPanel.add(chkUPGMA, gbc);
        gbc.gridy++;
        compressionPanel.add(chkNJ, gbc);
        gbc.gridx = 0;
        gbc.gridy++;
        compressionPanel.add(new JLabel("Group Kinship"), gbc);
        gbc.gridx++;
        compressionPanel.add(chkAvg, gbc);
        gbc.gridy++;
        compressionPanel.add(chkMin, gbc);
        gbc.gridy++;
        compressionPanel.add(chkMax, gbc);
        gbc.gridy++;
        compressionPanel.add(chkMedian, gbc);
        
        gbc.gridy = 0;
        gbc.gridx = 0;
        gbc.weightx = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(6,4,6,4);
        groupNumberPanel.add(new JLabel("Compression Levels"), gbc);
        gbc.gridx++;
        gbc.weightx = 0;
        gbc.insets = new Insets(6,4,6,4);
        groupNumberPanel.add(radioCompRange, gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.EAST;
        groupNumberPanel.add(new JLabel("From"), gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.WEST;
        groupNumberPanel.add(txtCompFrom, gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.EAST;
        groupNumberPanel.add(new JLabel("To"), gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.WEST;
        groupNumberPanel.add(txtCompTo, gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.EAST;
        groupNumberPanel.add(new JLabel("By"), gbc);
        gbc.gridx++;
        gbc.weightx = 1;
        gbc.anchor = GridBagConstraints.WEST;
        groupNumberPanel.add(txtCompBy, gbc);
        gbc.gridx = 1;
        gbc.gridy++;
        groupNumberPanel.add(radioCompList, gbc);
        gbc.gridx++;
        gbc.gridwidth = 6;
        groupNumberPanel.add(txtCompNumberList, gbc);

        gbc.gridy++;
        gbc.gridx = 0;
        gbc.gridwidth = 1;
        gbc.weightx = 0;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.insets = new Insets(6,4,6,4);
        groupNumberPanel.add(new JLabel("Group Sizes"), gbc);
        gbc.gridx++;
        gbc.weightx = 0;
        gbc.insets = new Insets(6,4,6,4);
        groupNumberPanel.add(radioGroupRange, gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.EAST;
        groupNumberPanel.add(new JLabel("From"), gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.WEST;
        groupNumberPanel.add(txtGroupFrom, gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.EAST;
        groupNumberPanel.add(new JLabel("To"), gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.WEST;
        groupNumberPanel.add(txtGroupTo, gbc);
        gbc.gridx++;
        gbc.anchor = GridBagConstraints.EAST;
        groupNumberPanel.add(new JLabel("By"), gbc);
        gbc.gridx++;
        gbc.weightx = 1;
        gbc.anchor = GridBagConstraints.WEST;
        groupNumberPanel.add(txtGroupBy, gbc);
        gbc.gridx = 1;
        gbc.gridy++;
        groupNumberPanel.add(radioGroupList, gbc);
        gbc.gridx++;
        gbc.gridwidth = 6;
        groupNumberPanel.add(txtGroupNumberList, gbc);

        
        centerPanel.add(compressionPanel, BorderLayout.NORTH);
        centerPanel.add(groupNumberPanel, BorderLayout.CENTER);
        theContentPane.add(centerPanel, BorderLayout.CENTER);
        
        JButton btnOK = new JButton("Run");
        btnOK.addActionListener(this);
        btnOK.setActionCommand("OK");
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(this);
        btnCancel.setActionCommand("Cancel");
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 0;
        gbc.gridwidth = 1;
        gbc.anchor = GridBagConstraints.CENTER;
        gbc.insets = new Insets(5,4,8,4);
        buttonPanel.add(btnOK, gbc);
        gbc.gridx++;
        buttonPanel.add(btnCancel, gbc);
        theContentPane.add(buttonPanel, BorderLayout.SOUTH);
        pack();
	}
    
	@Override
	public void actionPerformed(ActionEvent evt) {
		if (evt.getActionCommand().equals("OK")) {
			this.setVisible(false);
		} else if (evt.getActionCommand().equals("Cancel")) {
			this.setVisible(false);
		} else if (evt.getActionCommand().equals("nomarkers")) {
			if (chkNoMarkers.isSelected()) {
				chkP3D.setEnabled(false);
			} else {
				chkP3D.setEnabled(true);
			}
		} else if (evt.getActionCommand().equals("compress")) {
			if (chkCompression.isSelected()) {
				chkUPGMA.setEnabled(true);
				chkNJ.setEnabled(true);
				chkAvg.setEnabled(true);
				chkMin.setEnabled(true);
				chkMax.setEnabled(true);
				chkMedian.setEnabled(true);
				radioGroupRange.setEnabled(true);
				radioGroupList.setEnabled(true);
				txtGroupBy.setEnabled(true);
				txtGroupFrom.setEnabled(true);
				txtGroupNumberList.setEnabled(true);
				txtGroupTo.setEnabled(true);
//				radioComp.setEnabled(true);
//				txtCompressionFrom.setEnabled(true);
//				txtCompressionTo.setEnabled(true);
//				txtCompressionLevels.setEnabled(true);
			} else {
				chkUPGMA.setEnabled(false);
				chkNJ.setEnabled(false);
				chkAvg.setEnabled(false);
				chkMin.setEnabled(false);
				chkMax.setEnabled(false);
				chkMedian.setEnabled(false);
//				radioFromTo.setEnabled(false);
//				radioList.setEnabled(false);
//				txtBy.setEnabled(false);
//				txtFrom.setEnabled(false);
//				txtGroupNumberList.setEnabled(false);
//				txtTo.setEnabled(false);
//				radioComp.setEnabled(false);
//				txtCompressionFrom.setEnabled(false);
//				txtCompressionTo.setEnabled(false);
//				txtCompressionLevels.setEnabled(false);
			}
		}
		
	}
	
	ArrayList<Integer> getListOfGroups() {
		ArrayList<Integer> groupList = new ArrayList<Integer>();
		if (radioGroupList.isSelected()) {
			String[] groups = txtGroupNumberList.getText().split(",");
			int n = groups.length;
			try {
				for (int i = 0; i < n; i++) groupList.add(Integer.parseInt(groups[i]));
			} catch (Exception e) {
				String msg = "Illegal character in group list.";
				JOptionPane.showMessageDialog(getParent(), msg, "Illegal List", JOptionPane.ERROR_MESSAGE);
				return null;
			}
		}
		return null;
	}
	
	public static void main(String[] args) {
		MLMNewOptionDialog mod = new MLMNewOptionDialog(null);
		mod.setVisible(true);
		System.exit(-1);
	}
}

