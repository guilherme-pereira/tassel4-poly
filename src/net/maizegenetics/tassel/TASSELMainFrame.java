/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license, version 2 and without
 * any warranty or technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
//Title:      TASSELMainFrame
//Version:
//Copyright:  Copyright (c) 1997
//Author:     Ed Buckler
//Company:    NCSU
package net.maizegenetics.tassel;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.prefs.TasselPrefs;

import javax.swing.*;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.*;
import java.io.*;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import net.maizegenetics.baseplugins.ArchaeopteryxPlugin;
import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;
import net.maizegenetics.baseplugins.CreateTreePlugin;
import net.maizegenetics.baseplugins.ExportPlugin;
import net.maizegenetics.baseplugins.FileLoadPlugin;
import net.maizegenetics.baseplugins.FilterAlignmentPlugin;
import net.maizegenetics.baseplugins.FilterSiteNamePlugin;
import net.maizegenetics.baseplugins.FilterTaxaAlignmentPlugin;
import net.maizegenetics.baseplugins.FilterTaxaPropertiesPlugin;
import net.maizegenetics.baseplugins.FilterTraitsPlugin;
import net.maizegenetics.baseplugins.FixedEffectLMPlugin;
import net.maizegenetics.baseplugins.FlapjackLoadPlugin;
import net.maizegenetics.baseplugins.GenotypeImputationPlugin;
import net.maizegenetics.baseplugins.GenotypeSummaryPlugin;
import net.maizegenetics.baseplugins.Grid2dDisplayPlugin;
import net.maizegenetics.baseplugins.IntersectionAlignmentPlugin;
import net.maizegenetics.baseplugins.KinshipPlugin;
import net.maizegenetics.baseplugins.LinkageDiseqDisplayPlugin;
import net.maizegenetics.baseplugins.LinkageDisequilibriumPlugin;
import net.maizegenetics.baseplugins.MLMPlugin;
import net.maizegenetics.baseplugins.ManhattanDisplayPlugin;
import net.maizegenetics.baseplugins.MergeAlignmentsPlugin;
import net.maizegenetics.baseplugins.PlinkLoadPlugin;
import net.maizegenetics.baseplugins.QQDisplayPlugin;
import net.maizegenetics.baseplugins.SeparatePlugin;
import net.maizegenetics.baseplugins.SequenceDiversityPlugin;
import net.maizegenetics.baseplugins.SynonymizerPlugin;
import net.maizegenetics.baseplugins.TableDisplayPlugin;
import net.maizegenetics.baseplugins.TreeDisplayPlugin;
import net.maizegenetics.baseplugins.UnionAlignmentPlugin;
import net.maizegenetics.baseplugins.chart.ChartDisplayPlugin;
import net.maizegenetics.baseplugins.genomicselection.RidgeRegressionEmmaPlugin;
import net.maizegenetics.baseplugins.numericaltransform.NumericalTransformPlugin;
import net.maizegenetics.gui.PrintHeapAction;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.ThreadedPluginListener;
import net.maizegenetics.progress.ProgressPanel;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 * TASSELMainFrame
 *
 */
public class TASSELMainFrame extends JFrame implements ActionListener {

    private static final Logger myLogger = Logger.getLogger(TASSELMainFrame.class);
    public static final String version = "4.3.7";
    public static final String versionDate = "April 17, 2014";
    private DataTreePanel myDataTreePanel;
    private String tasselDataFile = "TasselDataFile";
    //a variable to control when the progress bar was last updated
    private JFileChooser filerSave = new JFileChooser();
    private JFileChooser filerOpen = new JFileChooser();
    private JScrollPane reportPanelScrollPane = new JScrollPane();
    private JTextArea reportPanelTextArea = new JTextArea();
    JScrollPane mainPanelScrollPane = new JScrollPane();
    JPanel mainDisplayPanel = new JPanel();
    private JTextArea mainPanelTextArea = new JTextArea();
    private JTextField myStatusTextField = new JTextField();
    private JMenuItem openCompleteDataTreeMenuItem = new JMenuItem();
    private JMenuItem openDataMenuItem = new JMenuItem();
    private JMenuItem saveCompleteDataTreeMenuItem = new JMenuItem();
    private JMenuItem saveDataTreeAsMenuItem = new JMenuItem();
    private PreferencesDialog thePreferencesDialog;
    private final ProgressPanel myProgressPanel = ProgressPanel.getInstance();
    private HashMap<JMenuItem, Plugin> myMenuItemHash = new HashMap<JMenuItem, Plugin>();

    public TASSELMainFrame() {
        try {
            loadSettings();
            myDataTreePanel = new DataTreePanel(this);
            myDataTreePanel.setToolTipText("Data Tree Panel");
            addMenuBar();
            initializeMyFrame();

            this.setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage) " + version);

            myLogger.info("Tassel Version: " + version + "  Date: " + versionDate);
            myLogger.info("Max Available Memory Reported by JVM: " + Utils.getMaxHeapSizeMB() + " MB");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    //Component initialization
    private void initializeMyFrame() throws Exception {
        getContentPane().setLayout(new BorderLayout());
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();

        // it is time for TASSEL to claim more (almost all) of the screen real estate for itself
        // this size was selected so as to encourage the user to resize to full screen,  thereby
        // insuring that all parts of the frame are visible.
        setSize(new Dimension(screenSize.width * 19 / 20, screenSize.height * 19 / 20));
        setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage)");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });

        filerSave.setDialogType(JFileChooser.SAVE_DIALOG);

        JSplitPane dataTreeReportPanelsSplitPanel = new JSplitPane();
        dataTreeReportPanelsSplitPanel.setOrientation(JSplitPane.VERTICAL_SPLIT);

        reportPanelTextArea.setEditable(false);
        reportPanelTextArea.setToolTipText("Report Panel");
        mainPanelTextArea.setDoubleBuffered(true);
        mainPanelTextArea.setEditable(false);
        mainPanelTextArea.setFont(new java.awt.Font("Monospaced", 0, 12));
        mainPanelTextArea.setToolTipText("Main Panel");

        myStatusTextField.setBackground(Color.lightGray);
        myStatusTextField.setBorder(null);
        saveCompleteDataTreeMenuItem.setText("Save Data Tree");
        saveCompleteDataTreeMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                saveCompleteDataTreeMenuItem_actionPerformed(e);
            }
        });
        saveDataTreeAsMenuItem.setText("Save Data Tree As ...");
        saveDataTreeAsMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                saveDataTreeMenuItem_actionPerformed(e);
            }
        });
        openCompleteDataTreeMenuItem.setText("Open Data Tree");
        openCompleteDataTreeMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {

                openCompleteDataTreeMenuItem_actionPerformed(e);
            }
        });

        openDataMenuItem.setText("Open Data Tree...");
        openDataMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {

                openDataMenuItem_actionPerformed(e);
            }
        });

        JSplitPane dataTreeReportMainPanelsSplitPanel = new JSplitPane();

        getContentPane().add(dataTreeReportMainPanelsSplitPanel, BorderLayout.CENTER);

        dataTreeReportMainPanelsSplitPanel.add(dataTreeReportPanelsSplitPanel, JSplitPane.LEFT);

        dataTreeReportPanelsSplitPanel.add(myDataTreePanel, JSplitPane.TOP);

        JSplitPane reportProgressSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

        JPanel reportPanel = new JPanel(new BorderLayout());

        reportProgressSplitPane.add(reportPanel, JSplitPane.TOP);

        reportPanel.add(reportPanelScrollPane, BorderLayout.CENTER);

        reportPanelScrollPane.getViewport().add(reportPanelTextArea, null);

        reportProgressSplitPane.add(new JScrollPane(myProgressPanel), JSplitPane.BOTTOM);

        dataTreeReportPanelsSplitPanel.add(reportProgressSplitPane, JSplitPane.BOTTOM);

        dataTreeReportMainPanelsSplitPanel.add(mainDisplayPanel, JSplitPane.RIGHT);

        mainDisplayPanel.setLayout(new BorderLayout());
        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);
        mainDisplayPanel.add(mainPanelScrollPane, BorderLayout.CENTER);

        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);

        getContentPane().add(myStatusTextField, BorderLayout.SOUTH);

        dataTreeReportMainPanelsSplitPanel.setDividerLocation(getSize().width / 4);

        dataTreeReportPanelsSplitPanel.setDividerLocation((int) (getSize().height / 3.5));

        reportProgressSplitPane.setDividerLocation((int) (getSize().height / 3.5));
    }

    private void addMenuBar() {

        JMenuBar jMenuBar = new JMenuBar();

        jMenuBar.add(getFileMenu());
        jMenuBar.add(getDataMenu());
        jMenuBar.add(getFiltersMenu());
        jMenuBar.add(getAnalysisMenu());
        jMenuBar.add(getResultsMenu());
        jMenuBar.add(getHelpMenu());

        this.setJMenuBar(jMenuBar);

    }

    //Help | About action performed
    private void helpAbout_actionPerformed(ActionEvent e) {
        AboutBox dlg = new AboutBox(this);
        Dimension dlgSize = dlg.getPreferredSize();
        Dimension frmSize = getSize();
        Point loc = getLocation();
        dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height - dlgSize.height) / 2 + loc.y);
        dlg.setModal(true);
        dlg.setVisible(true);
    }

    public void sendMessage(String text) {
        myStatusTextField.setForeground(Color.BLACK);
        myStatusTextField.setText(text);
    }

    public void sendErrorMessage(String text) {
        myStatusTextField.setForeground(Color.RED);
        myStatusTextField.setText(text);
    }

    public void setMainText(String text) {
        mainPanelTextArea.setText(text);
    }

    public void setNoteText(String text) {
        reportPanelTextArea.setText(text);
    }

    private void loadSettings() {
        filerOpen.setCurrentDirectory(new File(TasselPrefs.getOpenDir()));
        filerSave.setCurrentDirectory(new File(TasselPrefs.getSaveDir()));
    }

    /**
     * Provides a save filer that remembers the last location something was
     * saved to
     */
    private File getSaveFile() {

        File saveFile = null;
        int returnVal = filerSave.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {
            saveFile = filerSave.getSelectedFile();

            TasselPrefs.putSaveDir(filerSave.getCurrentDirectory().getPath());
        }
        return saveFile;
    }

    /**
     * Provides a open filer that remember the last location something was
     * opened from
     */
    private File getOpenFile() {

        File openFile = null;

        int returnVal = filerOpen.showOpenDialog(this);
        System.out.println("returnVal = " + returnVal);

        System.out.println("JFileChooser.OPEN_DIALOG " + JFileChooser.OPEN_DIALOG);

        if (returnVal == JFileChooser.OPEN_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            openFile = filerOpen.getSelectedFile();
            System.out.println("openFile = " + openFile);
            TasselPrefs.putOpenDir(filerOpen.getCurrentDirectory().getPath());
        }

        return openFile;
    }

    public void addDataSet(DataSet theDataSet, String defaultNode) {
        myDataTreePanel.addDataSet(theDataSet, defaultNode);
    }

    private void saveDataTree(String file) {

        Map dataToSerialize = new LinkedHashMap();
        Map dataFromTree = myDataTreePanel.getDataList();
        StringBuilder builder = new StringBuilder();

        Iterator itr = dataFromTree.keySet().iterator();
        while (itr.hasNext()) {

            Datum currentDatum = (Datum) itr.next();
            String currentNode = (String) dataFromTree.get(currentDatum);

            try {
                ByteArrayOutputStream out = new ByteArrayOutputStream();
                ObjectOutputStream oos = new ObjectOutputStream(out);
                oos.writeObject(currentDatum);
                oos.close();
                if (out.toByteArray().length > 0) {
                    dataToSerialize.put(currentDatum, currentNode);
                }
            } catch (Exception e) {
                myLogger.warn("saveDataTree: object: " + currentDatum.getName() + " type: " + currentDatum.getData().getClass().getName() + " does not serialize.");
                myLogger.warn("saveDataTree: message: " + e.getMessage());
                if (builder.length() == 0) {
                    builder.append("Due to error, these data sets could not\n");
                    builder.append("included in the saved file...");
                }
                builder.append("Data set: ");
                builder.append(currentDatum.getName());
                builder.append(" type: ");
                builder.append(currentDatum.getData().getClass().getName());
                builder.append("\n");
            }

        }

        try {

            File theFile = new File(Utils.addSuffixIfNeeded(file, ".zip"));
            FileOutputStream fos = new FileOutputStream(theFile);
            java.util.zip.ZipOutputStream zos = new ZipOutputStream(fos);

            ZipEntry thisEntry = new ZipEntry("DATA");
            zos.putNextEntry(thisEntry);
            ObjectOutputStream oos = new ObjectOutputStream(zos);
            oos.writeObject(dataToSerialize);
            oos.flush();
            zos.closeEntry();
            fos.close();
            sendMessage("Data saved to " + theFile.getAbsolutePath());

        } catch (Exception ee) {
            sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }

        if (builder.length() != 0) {
            JOptionPane.showMessageDialog(this, builder.toString(), "These data sets not saved...", JOptionPane.INFORMATION_MESSAGE);
        }

    }

    private boolean readDataTree(String file) {

        String dataTreeLoadFailed = "Unable to open the saved data tree.  The file format of this version is "
                + "incompatible with other versions.";

        boolean loadedDataTreePanel = false;
        try {

            FileInputStream fis = null;
            ObjectInputStream ois = null;
            if (file.endsWith("zip")) {

                fis = new FileInputStream(file);
                java.util.zip.ZipInputStream zis = new java.util.zip.ZipInputStream(fis);
                zis.getNextEntry();
                ois = new ObjectInputStream(zis);

            } else {

                fis = new FileInputStream(file);
                ois = new ObjectInputStream(fis);

            }

            try {
                Map data = (Map) ois.readObject();
                Iterator itr = data.keySet().iterator();
                while (itr.hasNext()) {
                    Datum currentDatum = (Datum) itr.next();
                    String currentNode = (String) data.get(currentDatum);
                    myDataTreePanel.addDatum(currentNode, currentDatum);
                }
                loadedDataTreePanel = true;
            } catch (InvalidClassException ice) {
                JOptionPane.showMessageDialog(this, dataTreeLoadFailed, "Incompatible File Format", JOptionPane.INFORMATION_MESSAGE);
            } finally {
                fis.close();
            }

            if (loadedDataTreePanel) {
                sendMessage("Data loaded.");
            }

        } catch (FileNotFoundException fnfe) {
            JOptionPane.showMessageDialog(this, "File not found: " + file, "File not found", JOptionPane.INFORMATION_MESSAGE);
            sendErrorMessage("Data tree could not be loaded.");
            return false;
        } catch (Exception ee) {
            JOptionPane.showMessageDialog(this, dataTreeLoadFailed + ee, "Incompatible File Format", JOptionPane.INFORMATION_MESSAGE);
            sendErrorMessage("Data tree could not be loaded.");
            return false;
        }

        return loadedDataTreePanel;

    }

    private void helpButton_actionPerformed(ActionEvent e) {
        HelpDialog theHelpDialog = new HelpDialog(this);
        theHelpDialog.setLocationRelativeTo(this);
        theHelpDialog.setVisible(true);
    }

    private void openCompleteDataTreeMenuItem_actionPerformed(ActionEvent e) {

        String dataFileName = tasselDataFile + ".zip";
        File dataFile = new File(dataFileName);
        if (dataFile.exists()) {
            readDataTree(dataFileName);
        } else if (new File("QPGADataFile").exists()) {
            // this exists to maintain backward compatibility with previous versions (pre-v0.99)
            readDataTree("QPGADataFile");
        } else {
            JOptionPane.showMessageDialog(this, "File: " + dataFile.getAbsolutePath() + " does not exist.\n"
                    + "Try using File/Open Data Tree...");
        }
    }

    private void openDataMenuItem_actionPerformed(ActionEvent e) {

        File f = getOpenFile();
        if (f != null) {
            readDataTree(f.getAbsolutePath());
        }
    }

    private void saveDataTreeMenuItem_actionPerformed(ActionEvent e) {

        File f = getSaveFile();
        if (f != null) {
            saveDataTree(f.getAbsolutePath());
        }
    }

    private void saveCompleteDataTreeMenuItem_actionPerformed(ActionEvent e) {
        saveDataTree(tasselDataFile + ".zip");
    }

    private void preferencesMenuItem_actionPerformed(ActionEvent e) {
        if (thePreferencesDialog == null) {
            thePreferencesDialog = new PreferencesDialog();
            thePreferencesDialog.pack();
        }

        thePreferencesDialog.setLocationRelativeTo(this);
        thePreferencesDialog.setVisible(true);
    }

    public void updateMainDisplayPanel(JPanel panel) {
        mainDisplayPanel.removeAll();
        mainDisplayPanel.add(panel, BorderLayout.CENTER);
        mainDisplayPanel.repaint();
        mainDisplayPanel.validate();
    }

    public DataTreePanel getDataTreePanel() {
        return myDataTreePanel;
    }

    public ProgressPanel getProgressPanel() {
        return myProgressPanel;
    }

    private JMenuItem createMenuItem(Plugin theTP) {
        return createMenuItem(theTP, -1);
    }

    private JMenuItem createMenuItem(Plugin theTP, int mnemonic) {
        ImageIcon icon = theTP.getIcon();
        JMenuItem menuItem = new JMenuItem(theTP.getButtonName(), icon);
        if (mnemonic != -1) {
            menuItem.setMnemonic(mnemonic);
        }
        int pixels = 30;
        if (icon != null) {
            pixels -= icon.getIconWidth();
            pixels /= 2;
        }
        menuItem.setIconTextGap(pixels);
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.setToolTipText(theTP.getToolTipText());
        menuItem.addActionListener(this);
        theTP.addListener(myDataTreePanel);
        myMenuItemHash.put(menuItem, theTP);
        return menuItem;
    }

    private JMenu getFiltersMenu() {
        JMenu result = new JMenu("Filter");
        result.setMnemonic(KeyEvent.VK_F);
        result.add(createMenuItem(new FilterAlignmentPlugin(this, true)));
        result.add(createMenuItem(new FilterSiteNamePlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaAlignmentPlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaPropertiesPlugin(this, true)));
        result.add(createMenuItem(new FilterTraitsPlugin(this, true)));
        return result;
    }

    private JMenu getDataMenu() {

        JMenu result = new JMenu("Data");
        result.setMnemonic(KeyEvent.VK_D);

        PlinkLoadPlugin plinkLoadPlugin = new PlinkLoadPlugin(this, true);
        plinkLoadPlugin.addListener(myDataTreePanel);

        FlapjackLoadPlugin flapjackLoadPlugin = new FlapjackLoadPlugin(this, true);
        flapjackLoadPlugin.addListener(myDataTreePanel);

        result.add(createMenuItem(new FileLoadPlugin(this, true, plinkLoadPlugin, flapjackLoadPlugin), KeyEvent.VK_L));
        result.add(createMenuItem(new ExportPlugin(this, true)));
        result.add(createMenuItem(new ConvertSBitTBitPlugin(this, true)));
        result.add(createMenuItem(new GenotypeImputationPlugin(this, true)));
        result.add(createMenuItem(new NumericalTransformPlugin(this, true)));
        result.add(createMenuItem(new SynonymizerPlugin(this, true)));
        result.add(createMenuItem(new IntersectionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new UnionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new MergeAlignmentsPlugin(this, true)));
        result.add(createMenuItem(new SeparatePlugin(this, true)));
        result.addSeparator();

        JMenuItem delete = new JMenuItem("Delete Dataset");
        delete.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                myDataTreePanel.deleteSelectedNodes();
            }
        });
        delete.setToolTipText("Delete Dataset");
        URL imageURL = TASSELMainFrame.class.getResource("images/trash.gif");
        if (imageURL == null) {
            delete.setIconTextGap(30);
        } else {
            delete.setIcon(new ImageIcon(imageURL));
            delete.setIconTextGap(6);
        }
        delete.setIconTextGap(6);
        result.add(delete);

        return result;
    }

    private JMenu getAnalysisMenu() {

        JMenu result = new JMenu("Analysis");
        result.setMnemonic(KeyEvent.VK_A);

        result.add(createMenuItem(new SequenceDiversityPlugin(this, true)));
        result.add(createMenuItem(new LinkageDisequilibriumPlugin(this, true)));
        result.add(createMenuItem(new CreateTreePlugin(this, true)));
        result.add(createMenuItem(new KinshipPlugin(this, true)));
        result.add(createMenuItem(new FixedEffectLMPlugin(this, true)));
        result.add(createMenuItem(new MLMPlugin(this, true)));
        result.add(createMenuItem(new RidgeRegressionEmmaPlugin(this, true)));
        result.add(createMenuItem(new GenotypeSummaryPlugin(this, true)));
//        result.add(createMenuItem(new StepwiseOLSModelFitterPlugin(this, true)));
        return result;
    }

    private JMenu getResultsMenu() {

        JMenu result = new JMenu("Results");
        result.setMnemonic(KeyEvent.VK_R);

        result.add(createMenuItem(new TableDisplayPlugin(this, true)));
        result.add(createMenuItem(new ArchaeopteryxPlugin(this, true)));
        result.add(createMenuItem(new Grid2dDisplayPlugin(this, true)));
        result.add(createMenuItem(new LinkageDiseqDisplayPlugin(this, true)));
        result.add(createMenuItem(new ChartDisplayPlugin(this, true)));
        result.add(createMenuItem(new QQDisplayPlugin(this, true)));
        result.add(createMenuItem(new ManhattanDisplayPlugin(this, true)));
        result.add(createMenuItem(new TreeDisplayPlugin(this, true)));
        return result;

    }

    private JMenu getFileMenu() {
        JMenu fileMenu = new JMenu();
        fileMenu.setText("File");
        fileMenu.add(saveCompleteDataTreeMenuItem);
        fileMenu.add(openCompleteDataTreeMenuItem);
        fileMenu.add(saveDataTreeAsMenuItem);
        fileMenu.add(openDataMenuItem);
        JMenuItem preferencesMenuItem = new JMenuItem();
        preferencesMenuItem.setText("Set Preferences");

        preferencesMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {

                preferencesMenuItem_actionPerformed(e);

            }
        });
        fileMenu.add(preferencesMenuItem);
        fileMenu.addSeparator();

        JMenuItem exitMenuItem = new JMenuItem("Exit");
        exitMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.exit(0);
            }
        });
        fileMenu.add(exitMenuItem);

        return fileMenu;
    }

    private JMenu getHelpMenu() {
        JMenu helpMenu = new JMenu();
        helpMenu.setMnemonic(KeyEvent.VK_H);
        helpMenu.setText("Help");

        JMenuItem helpMenuItem = new JMenuItem("Help Manual");
        helpMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                helpButton_actionPerformed(e);
            }
        });
        helpMenu.add(helpMenuItem);

        JMenuItem aboutMenuItem = new JMenuItem("About");
        aboutMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                helpAbout_actionPerformed(e);
            }
        });
        helpMenu.add(aboutMenuItem);

        JMenuItem memoryUsage = new JMenuItem("Show Memory");
        memoryUsage.addActionListener(PrintHeapAction.getInstance(this));
        memoryUsage.setToolTipText("Show Memory Usage");
        helpMenu.add(memoryUsage);

        return helpMenu;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        JMenuItem theMenuItem = (JMenuItem) e.getSource();
        Plugin theTP = this.myMenuItemHash.get(theMenuItem);
        PluginEvent event = new PluginEvent(myDataTreePanel.getSelectedTasselDataSet());
        ProgressPanel progressPanel = getProgressPanel();
        progressPanel.addPlugin(theTP);
        ThreadedPluginListener thread = new ThreadedPluginListener(theTP, event);
        thread.start();
    }
}
