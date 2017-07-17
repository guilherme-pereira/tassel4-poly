/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
package net.maizegenetics.tassel;

import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.ids.IdentifierSynonymizer;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.tree.Tree;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginListener;

import net.maizegenetics.baseplugins.FixedEffectLMPlugin;
import net.maizegenetics.baseplugins.SequenceDiversityPlugin;
import net.maizegenetics.baseplugins.LinkageDisequilibriumPlugin;
import net.maizegenetics.baseplugins.MLMPlugin;

import net.maizegenetics.util.Utils;

import javax.swing.*;
import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.*;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

import java.io.Serializable;

import java.net.URL;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.batik.util.gui.MemoryMonitor;

import org.apache.log4j.Logger;

public class DataTreePanel extends JPanel implements PluginListener {

    private static final Logger myLogger = Logger.getLogger(DataTreePanel.class);
    public static final String NODE_TYPE_DATA = "Data";
    public static final String NODE_TYPE_RESULT = "Result";
    public static final String NODE_TYPE_SEQUENCE = "Sequence";
    public static final String NODE_TYPE_POLYMORPHISMS = "Polymorphisms";
    public static final String NODE_TYPE_NUMERICAL = "Numerical";
    public static final String NODE_TYPE_MATRIX = "Matrix";
    public static final String NODE_TYPE_TREE = "Tree";
    public static final String NODE_TYPE_FUSIONS = "Fusions";
    public static final String NODE_TYPE_SYNONYMS = "Synonyms";
    public static final String NODE_TYPE_DIVERSITY = "Diversity";
    public static final String NODE_TYPE_SNP_ASSAYS = "SNP Assays";
    public static final String NODE_TYPE_LD = "LD";
    public static final String NODE_TYPE_ASSOCIATIONS = "Association";
    public static final String NODE_TYPE_VARIANCES = "Variances";
    public static final String NODE_TYPE_SYNONYMIZER = "Synonymizer";
    public static final String NODE_TYPE_STEPWISE = "Stepwise";
    public static final String NODE_TYPE_DEFAULT = NODE_TYPE_DATA;
    //Possible line styles...
    //"Angled", "Horizontal", and "None" (the default).
    private String myLineStyle = "Angled";
    private JTree myTree;
    private DefaultTreeModel myTreeModel;
    private TASSELMainFrame myTASSELMainFrame;
    private HashMap myNodeHash = new HashMap();
    private LinkedHashMap myDataSetList = new LinkedHashMap();
    private Datum myLastBookSelected;
    private DefaultMutableTreeNode myTopNode;
    private DefaultMutableTreeNode myDataNode;
    private DefaultMutableTreeNode myResultNode;
    private DefaultMutableTreeNode myLastNode = null;
    private TreeSelectionListener myTreeSelectionListener = null;

    public DataTreePanel(TASSELMainFrame theQAF) {
        super();

        myTASSELMainFrame = theQAF;

        myTopNode = new DefaultMutableTreeNode("top");
        createNodes();

        myTreeModel = new DefaultTreeModel(myTopNode);

        myTree = new JTree(myTreeModel);
        myTree.setRootVisible(false);
        myTree.setEditable(true);
        myTreeModel.addTreeModelListener(new MyTreeModelListener(this));
        myTree.getSelectionModel().setSelectionMode(TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION);
        myTree.setLargeModel(true);
        //Listen for when the selection changes.
        initSelectionListener();
        initKeyStrokeListener();

        myTree.putClientProperty("JTree.lineStyle", myLineStyle);

        URL tsBitURL = DataTreePanel.class.getResource("images/tsBit.gif");
        final ImageIcon tsBitIcon;
        if (tsBitURL != null) {
            tsBitIcon = new ImageIcon(tsBitURL);
        } else {
            tsBitIcon = null;
        }

        URL sBitURL = DataTreePanel.class.getResource("images/sBit.gif");
        final ImageIcon sBitIcon;
        if (sBitURL != null) {
            sBitIcon = new ImageIcon(sBitURL);
        } else {
            sBitIcon = null;
        }

        URL tBitURL = DataTreePanel.class.getResource("images/tBit.gif");
        final ImageIcon tBitIcon;
        if (tBitURL != null) {
            tBitIcon = new ImageIcon(tBitURL);
        } else {
            tBitIcon = null;
        }

        URL combineURL = DataTreePanel.class.getResource("images/combineAlign.gif");
        final ImageIcon combineIcon;
        if (combineURL != null) {
            combineIcon = new ImageIcon(combineURL);
        } else {
            combineIcon = null;
        }

        URL filterURL = DataTreePanel.class.getResource("images/filterAlign.gif");
        final ImageIcon filterIcon;
        if (filterURL != null) {
            filterIcon = new ImageIcon(filterURL);
        } else {
            filterIcon = null;
        }

        URL sBitCombineURL = DataTreePanel.class.getResource("images/sBitCombine.gif");
        final ImageIcon sBitCombineIcon;
        if (sBitCombineURL != null) {
            sBitCombineIcon = new ImageIcon(sBitCombineURL);
        } else {
            sBitCombineIcon = null;
        }

        URL sBitFilterURL = DataTreePanel.class.getResource("images/sBitFilter.gif");
        final ImageIcon sBitFilterIcon;
        if (sBitFilterURL != null) {
            sBitFilterIcon = new ImageIcon(sBitFilterURL);
        } else {
            sBitFilterIcon = null;
        }

        URL tBitFilterURL = DataTreePanel.class.getResource("images/tBitFilter.gif");
        final ImageIcon tBitFilterIcon;
        if (tBitFilterURL != null) {
            tBitFilterIcon = new ImageIcon(tBitFilterURL);
        } else {
            tBitFilterIcon = null;
        }


        myTree.setCellRenderer(new DefaultTreeCellRenderer() {
            public Component getTreeCellRendererComponent(JTree pTree,
                    Object pValue, boolean pIsSelected, boolean pIsExpanded,
                    boolean pIsLeaf, int pRow, boolean pHasFocus) {
                DefaultMutableTreeNode node = (DefaultMutableTreeNode) pValue;
                Datum nodeInfo = (Datum) node.getUserObject();
                Component result = super.getTreeCellRendererComponent(pTree, pValue, pIsSelected,
                        pIsExpanded, pIsLeaf, pRow, pHasFocus);
                Object data = nodeInfo.getData();
                if (data instanceof AlignmentMask) {
                    result.setForeground(((AlignmentMask) nodeInfo.getData()).getColor());
                }

                if (data instanceof Alignment) {
                    Alignment align = (Alignment) data;
                    if ((align.isTBitFriendly()) && (align.isSBitFriendly())) {
                        setIcon(tsBitIcon);
                    } else if (align.isSBitFriendly()) {
                        if ((align instanceof FilterAlignment) && (sBitFilterIcon != null)) {
                            setIcon(sBitFilterIcon);
                        } else if ((align instanceof CombineAlignment) && (sBitCombineIcon != null)) {
                            setIcon(sBitCombineIcon);
                        } else if (sBitIcon != null) {
                            setIcon(sBitIcon);
                        }
                    } else if (align.isTBitFriendly()) {
                        if ((align instanceof FilterAlignment) && (tBitFilterIcon != null)) {
                            setIcon(tBitFilterIcon);
                        } else if (tBitIcon != null) {
                            setIcon(tBitIcon);
                        }
                    } else if (align instanceof FilterAlignment) {
                        if (filterIcon != null) {
                            setIcon(filterIcon);
                        }
                    } else if (align instanceof CombineAlignment) {
                        if (combineIcon != null) {
                            setIcon(combineIcon);
                        }
                    }
                }
                return result;
            }
        });


        try {
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public DataTreePanel getThis() {
        return this;
    }

    private void processMyKeyStroke(KeyEvent e) {

        // only respond if the KeyEvent is a key released event.
        if (e.getID() == KeyEvent.KEY_RELEASED) {
            if (KeyEvent.VK_DELETE == e.getKeyCode()) {

                TreePath[] currentSelection = myTree.getSelectionPaths();
                String pluralQuestion = "Are you sure you want to delete these " + currentSelection.length + " nodes?";
                String question = "Are you sure you want to delete this node?";

                if (currentSelection.length > 1) {
                    question = pluralQuestion;
                }
                int userResponse = JOptionPane.showOptionDialog(myTASSELMainFrame, question, "Node Deletion",
                        JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);

                if (userResponse == JOptionPane.YES_OPTION) {
                    deleteSelectedNodes();
                }
            }
            if (KeyEvent.VK_F12 == e.getKeyCode()) {    // open a memory monitor to show usage
                MemoryMonitor memMon = new MemoryMonitor();
                memMon.pack();
                memMon.setSize(new Dimension(250, 300));
                memMon.setVisible(true);
            }
        }
    }

    private void initKeyStrokeListener() {

        myTree.addKeyListener(new KeyAdapter() {
            /**
             * Invoked when a key has been released.
             */
            public void keyReleased(KeyEvent e) {
                super.keyReleased(e);
                processMyKeyStroke(e);
            }
        });
    }

    private void initSelectionListener() {

        myTreeSelectionListener = (new TreeSelectionListener() {
            public void valueChanged(TreeSelectionEvent e) {
                try {
                    DefaultMutableTreeNode node = null;
                    if ((e != null) && (e.getSource() instanceof DefaultMutableTreeNode)) {
                        node = (DefaultMutableTreeNode) e.getSource();
                    } else {
                        node = (DefaultMutableTreeNode) myTree.getLastSelectedPathComponent();
                    }

                    if (node == null) {
                        return;
                    }
                    if (myLastNode == node) {
                        return;
                    }
                    myLastNode = node;
                    Object nodeInfo = node.getUserObject();
                    myLastBookSelected = (Datum) nodeInfo;

                    Datum book = (Datum) nodeInfo;
                    myLogger.info("initSelectionListener: node type: " + book.getDataType());
                    if (!(book.getData() instanceof Serializable)) {
                        myLogger.info("This is not serializable.");
                    }
                    StringBuilder builder = new StringBuilder();
                    if (book.getData() instanceof Alignment) {
                        Alignment a = (Alignment) book.getData();
                        builder.append("Number of sequences: ");
                        builder.append(a.getSequenceCount());
                        builder.append("\n");
                        builder.append("Number of sites: ");
                        builder.append(a.getSiteCount());
                        builder.append("\n");
                        Locus[] loci = a.getLoci();
                        boolean first = true;
                        int numLoci = 0;
                        if (loci != null) {
                            numLoci = loci.length;
                        }
                        for (int i = 0; i < numLoci; i++) {
                            String name = loci[i].getName();
                            if ((name == null) || (name.length() == 0)) {
                                continue;
                            }
                            if (first) {
                                builder.append("Loci: ");
                                first = false;
                            } else {
                                builder.append(", ");
                            }
                            builder.append(name);
                        }
                        if (!first) {
                            builder.append("\n");
                        }
                    }
                    if (book.getData() instanceof FilterAlignment) {
                        FilterAlignment a = (FilterAlignment) book.getData();
                        builder.append("FilterAlignment...\n");
                        if (a.isSiteFilter()) {
                            builder.append("Site Filter\n");
                        }
                        if (a.isSiteFilterByRange()) {
                            builder.append("Site Range Filter\n");
                        }
                        if (a.isTaxaFilter()) {
                            builder.append("Taxa Filter\n");
                        }
                        builder.append("Base Type: ");
                        builder.append(Utils.getBasename(a.getBaseAlignment().getClass().getName()));
                        builder.append("\n");
                    }
                    if (book.getData() instanceof TableReport) {
                        TableReport tr = (TableReport) book.getData();
                        builder.append("Table Title: ");
                        builder.append(tr.getTableTitle());
                        builder.append("\n");
                        builder.append("Number of columns: ");
                        builder.append(tr.getColumnCount());
                        builder.append("\n");
                        builder.append("Number of rows: ");
                        builder.append(tr.getRowCount());
                        builder.append("\n");
                        builder.append("Number of elements: ");
                        builder.append(tr.getElementCount());
                        builder.append("\n");
                    }
                    if (book.getData() instanceof AlignmentMask) {
                        AlignmentMask mask = (AlignmentMask) book.getData();
                    }

                    String comment = book.getComment();
                    if ((comment != null) && (comment.length() != 0)) {
                        builder.append(comment);
                    }
                    myTASSELMainFrame.setNoteText(builder.toString());
                    if (book.getData() != null) {
                        myTASSELMainFrame.mainDisplayPanel.removeAll();

                        if (book.getData() instanceof TableReport) {
                            //This method issues that giant files are not opened as JTables
                            //In the future it may be good to add a getSize method to TableReport
                            int size = ((TableReport) book.getData()).getElementCount();
                            myLogger.info("initSelectionListener: Table Report Size: " + size);
                            if (size == 0) {
                                JPanel blankPanel = new JPanel();
                                blankPanel.setLayout(new BorderLayout());
                                blankPanel.add(new JLabel("     Nothing to Display"), BorderLayout.CENTER);
                                myTASSELMainFrame.mainDisplayPanel.add(blankPanel, BorderLayout.CENTER);
                            } else {
                                TableReportPanel theATP;
                                if (book.getData() instanceof Phenotype) {
                                    theATP = new TableReportPanel((Phenotype) book.getData());
                                } else {
                                    theATP = new TableReportPanel((TableReport) book.getData());
                                }
                                myTASSELMainFrame.mainDisplayPanel.add(theATP, BorderLayout.CENTER);
                            }
                        } else if (book.getData() instanceof Alignment) {
                            Alignment align = (Alignment) book.getData();
                            List masks = new ArrayList();
                            for (int i = 0, n = node.getChildCount(); i < n; i++) {
                                try {
                                    DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) node.getChildAt(i);
                                    Datum currentDatum = (Datum) currentNode.getUserObject();
                                    Object currentMask = currentDatum.getData();
                                    if ((currentMask != null) && (currentMask instanceof AlignmentMask)) {
                                        masks.add(currentMask);
                                    }
                                } catch (Exception ex) {
                                    ex.printStackTrace();
                                }
                            }

                            SeqViewerPanel seqViewer = null;
                            if (masks.size() != 0) {
                                AlignmentMask[] masksArray = new AlignmentMask[masks.size()];
                                masks.toArray(masksArray);
                                seqViewer = SeqViewerPanel.getInstance(align, masksArray, getThis());
                            } else {
                                seqViewer = SeqViewerPanel.getInstance(align, getThis());
                            }
                            myTASSELMainFrame.mainDisplayPanel.add(seqViewer, BorderLayout.CENTER);
                        } else if (book.getData() instanceof AlignmentMask) {
                            AlignmentMask mask = (AlignmentMask) book.getData();
                            if (mask.getColor() != null) {
                                Color tempColor = JColorChooser.showDialog(myTASSELMainFrame, "Select Mask Color...", mask.getColor());
                                //myTree.setSelectionPath(new TreePath(parentNode.getPath()));
                                if (tempColor != null) {
                                    mask.setColor(tempColor);
                                    book.setName(mask.toString());
                                }
                            }
                            myLastNode = null;
                            DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) node.getParent();
                            TreeSelectionEvent event = new TreeSelectionEvent(parentNode, null, false, null, null);
                            myTreeSelectionListener.valueChanged(event);
                        } else {
                            myTASSELMainFrame.mainDisplayPanel.add(myTASSELMainFrame.mainPanelScrollPane, BorderLayout.CENTER);
                            //x = book.getData().toString();
                            String s = book.getData().toString();
                            if (s.length() > 2000000) {
                                s = s.substring(1, 2000000) + "\n Truncated view.  Too much to display.  Save it to a file.";
                            }
                            myTASSELMainFrame.setMainText(s);
                        }

                        myTASSELMainFrame.sendMessage(book.getData().getClass().toString());

                        myTASSELMainFrame.mainDisplayPanel.revalidate();
                        myTASSELMainFrame.mainDisplayPanel.repaint();

                    }

                } catch (Throwable e1) {

                    String userMessage = "TASSEL has experienced an error.  ";
                    if (e1 instanceof java.lang.OutOfMemoryError) {
                        userMessage = "You have used up all of the memory allocated to the Java Virtual Machine.  \n"
                                + "It is recommneded that you adjust your heap settings and possibly add more memory to the computer\n"
                                + "Additionally, some operations are not recommended on a full dataset, i.e., select SNPs *before* determining LD";
                    }
                    JOptionPane.showMessageDialog(myTASSELMainFrame, userMessage, "Fatal Error", JOptionPane.ERROR_MESSAGE);
                    e1.printStackTrace();
                }
            }
        });

        myTree.addTreeSelectionListener(myTreeSelectionListener);

    }

    private void createNodes() {
        DefaultMutableTreeNode book = null;

        myDataNode = new DefaultMutableTreeNode(new Datum(NODE_TYPE_DATA, "Node on data tree", "Holds the basic data structures"));
        myTopNode.add(myDataNode);
        myNodeHash.put(NODE_TYPE_DATA, myDataNode);
        myResultNode = new DefaultMutableTreeNode(new Datum(NODE_TYPE_RESULT, "Node on data tree", "Holds the basic result structures"));
        myTopNode.add(myResultNode);
        myNodeHash.put(NODE_TYPE_RESULT, myResultNode);
        DefaultMutableTreeNode geneNode = new DefaultMutableTreeNode(new Datum(NODE_TYPE_SEQUENCE, "Node on data tree", "Please load some genes"));
        myDataNode.add(geneNode);
        myNodeHash.put(NODE_TYPE_SEQUENCE, geneNode);
        DefaultMutableTreeNode polymorphismNode = new DefaultMutableTreeNode(new Datum(NODE_TYPE_POLYMORPHISMS, "Node on data tree", "Please load some polymorphisms"));
        myDataNode.add(polymorphismNode);
        myNodeHash.put(NODE_TYPE_POLYMORPHISMS, polymorphismNode);
        DefaultMutableTreeNode phenotypeNode = new DefaultMutableTreeNode(new Datum(NODE_TYPE_NUMERICAL, "Node on data tree", "Please load some phenotypes"));
        myDataNode.add(phenotypeNode);
        myNodeHash.put(NODE_TYPE_NUMERICAL, phenotypeNode);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_MATRIX, "Node on data tree", "Kinship matrix"));
        myDataNode.add(book);
        myNodeHash.put(NODE_TYPE_MATRIX, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_TREE, "Node on data tree", "Cladograms and Trees"));
        myDataNode.add(book);
        myNodeHash.put(NODE_TYPE_TREE, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_FUSIONS, "Node on data tree", "Fusions between genes and datatypes"));
        myDataNode.add(book);
        myNodeHash.put(NODE_TYPE_FUSIONS, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_SYNONYMIZER, "Node on data tree", "Taxa Synonyms"));
        myDataNode.add(book);
        myNodeHash.put(NODE_TYPE_SYNONYMIZER, book);

        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_DIVERSITY, "Node on data tree", "Diversity"));
        myResultNode.add(book);
        myNodeHash.put(NODE_TYPE_DIVERSITY, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_SNP_ASSAYS, "Node on data tree", "SNP Extracted Data"));
        myResultNode.add(book);
        myNodeHash.put(NODE_TYPE_SNP_ASSAYS, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_LD, "Node on data tree", "Linkage Disequilibrium"));
        myResultNode.add(book);
        myNodeHash.put(NODE_TYPE_LD, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_ASSOCIATIONS, "Node on data tree", "Phenotypic Associations"));
        myResultNode.add(book);
        myNodeHash.put(NODE_TYPE_ASSOCIATIONS, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_VARIANCES, "Node on data tree", "Additive genetic and residual variance"));
        myResultNode.add(book);
        myNodeHash.put(NODE_TYPE_VARIANCES, book);
        book = new DefaultMutableTreeNode(new Datum(NODE_TYPE_STEPWISE, "Node on data tree", "Stepwise Regression"));
        myResultNode.add(book);
        myNodeHash.put(NODE_TYPE_STEPWISE, book);

    }

    private void jbInit() throws Exception {
        JScrollPane jScrollPane1 = new JScrollPane(myTree);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        setLayout(new BorderLayout());
        add(jScrollPane1, BorderLayout.CENTER);
    }

    private void deleteDatumFromList(Datum datum) {

        Iterator itr = myDataSetList.keySet().iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();
            if (current == datum) {
                itr.remove();
                break;
            }
        }

    }

    public void addDataSet(DataSet theDataSet, String defaultNode) {

        Plugin theCreator = theDataSet.getCreator();
        for (int i = 0; i < theDataSet.getSize(); i++) {
            Datum d = theDataSet.getData(i);
            if ((theCreator instanceof MLMPlugin)
                    || (theCreator instanceof FixedEffectLMPlugin)) {
                addDatum(NODE_TYPE_ASSOCIATIONS, d);
                continue;

            }


            if (theCreator instanceof SequenceDiversityPlugin) {
                addDatum(NODE_TYPE_DIVERSITY, d);
                continue;

            }


            if (theCreator instanceof LinkageDisequilibriumPlugin) {
                addDatum(NODE_TYPE_LD, d);
                continue;

            }


            if (d.getData() instanceof Alignment) {
                addDatum(NODE_TYPE_SEQUENCE, d);
                continue;

            }


            if (d.getData() instanceof AlignmentMask) {
                addDatum(d);
                continue;

            }


            if (d.getData() instanceof IdentifierSynonymizer) {
                addDatum(NODE_TYPE_SYNONYMIZER, d);
                continue;

            }


            if (d.getData() instanceof Phenotype) {
                addDatum(NODE_TYPE_NUMERICAL, d);
                continue;

            }


            if (d.getData() instanceof DistanceMatrix) {
                addDatum(NODE_TYPE_MATRIX, d);
                continue;

            }


            if (d.getData() instanceof TableReport) {
                addDatum(NODE_TYPE_NUMERICAL, d);
                continue;

            }


            if (d.getData() instanceof Tree) {
                addDatum(NODE_TYPE_TREE, d);
                continue;

            }


            if (defaultNode == null) {
                addDatum(NODE_TYPE_DEFAULT, d);
                continue;

            }



            addDatum(defaultNode, d);
        }

    }

    public synchronized void addDatum(String dataParent, Datum theDatum) {
        if (theDatum.getData() instanceof AlignmentMask) {
            addDatum(theDatum);
            return;
        }

        DefaultMutableTreeNode temp = findOnTree(theDatum);
        if (temp != null) {
            myTree.scrollPathToVisible(new TreePath(temp.getPath()));
            myTree.validate();
            myTree.repaint();
            return;
        }

        DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) myNodeHash.get(dataParent);
        DefaultMutableTreeNode childNode = new DefaultMutableTreeNode(theDatum);
        myNodeHash.put(theDatum.getName(), childNode);
        myTreeModel.insertNodeInto(childNode, parentNode, parentNode.getChildCount());
        myTree.scrollPathToVisible(new TreePath(childNode.getPath()));
        if (theDatum instanceof Serializable) {
            myDataSetList.put(theDatum, dataParent);
        }
    }

    public synchronized void addDatum(Datum theDatum) {
        AlignmentMask mask = null;
        try {
            mask = (AlignmentMask) theDatum.getData();
        } catch (Exception e) {
            throw new IllegalArgumentException("DataTreePanel: addDatum: input must be AlignmentMask.");
        }
        Alignment align = mask.getAlignment();
        Iterator itr = myNodeHash.values().iterator();
        while (itr.hasNext()) {
            DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) itr.next();
            Object current = ((Datum) parentNode.getUserObject()).getData();
            if ((current instanceof Alignment) && (current == align)) {
                DefaultMutableTreeNode childNode = new DefaultMutableTreeNode(theDatum);
                myTreeModel.insertNodeInto(childNode, parentNode, parentNode.getChildCount());
                myTree.scrollPathToVisible(new TreePath(childNode.getPath()));
                myLastNode = null;
                myTree.setSelectionPath(new TreePath(parentNode.getPath()));
                myTree.scrollPathToVisible(new TreePath(parentNode.getPath()));
                myTreeSelectionListener.valueChanged(null);
                if (theDatum instanceof Serializable) {
                    myDataSetList.put(theDatum, "NA");
                }
                break;
            }
        }
    }

    private DefaultMutableTreeNode findOnTree(Datum theDatum) {

        Object obj = theDatum.getData();
        Iterator itr = myNodeHash.values().iterator();
        while (itr.hasNext()) {
            DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) itr.next();
            Object current = ((Datum) parentNode.getUserObject()).getData();
            if (current == obj) {
                return parentNode;
            }
        }
        return null;

    }

    public DataSet getSelectedTasselDataSet() {
        int n = myTree.getSelectionCount();
        String parentNodeName = null;
        ArrayList<Datum> data = new ArrayList<Datum>();
        TreePath[] tp = myTree.getSelectionPaths();
        for (int i = 0; i < n; i++) {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) tp[i].getLastPathComponent();
            if (node == null) {
                return null;
            }
            //if the parents don't share the same name then leave blank

            if (i == 0) {
                parentNodeName = node.getParent().toString();
            } else if (!node.getParent().toString().equals(parentNodeName)) {
                parentNodeName = null;
            }

            data.add((Datum) node.getUserObject());
        }

        return new DataSet(data, null);
    }

    public void deleteSelectedNodes() {

        TreePath[] currentSelection = myTree.getSelectionPaths();
        if (currentSelection != null) {
            for (int i = 0; i < currentSelection.length; i++) {
                myLogger.info("Start Deleting at Selection: " + currentSelection[i]);
                if (currentSelection[i] != null) {
                    DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) (currentSelection[i].getLastPathComponent());
                    myLogger.info("Deleting Node : " + currentNode);
                    deleteNode(currentNode);
                }

            }
        }

    }

    private void deleteNode(DefaultMutableTreeNode currentNode) {

        int size = currentNode.getChildCount();
        for (int i = size - 1; i >= 0; i--) {
            DefaultMutableTreeNode current = (DefaultMutableTreeNode) currentNode.getChildAt(i);
            deleteNode(current);
        }

        if ((currentNode == myDataNode) || (currentNode == myResultNode)) {
            return;
        }

        MutableTreeNode parent = (MutableTreeNode) (currentNode.getParent());
        if ((parent != null) && (parent != myDataNode) && (parent != myResultNode)) {
            myTreeModel.removeNodeFromParent(currentNode);
            Iterator itr = myNodeHash.keySet().iterator();
            while (itr.hasNext()) {
                Object key = itr.next();
                if (myNodeHash.get(key) == currentNode) {
                    itr.remove();
                }

            }

            try {
                Datum datum = (Datum) currentNode.getUserObject();
                deleteDatumFromList(datum);
            } catch (Exception ex) {
                // do nothing
            }

            try {
                if (((Datum) currentNode.getUserObject()).getData() instanceof AlignmentMask) {
                    myLastNode = null;
                    myTree.setSelectionPath(new TreePath(((DefaultMutableTreeNode) parent).getPath()));
                    myTree.scrollPathToVisible(new TreePath(((DefaultMutableTreeNode) parent).getPath()));
                    myTreeSelectionListener.valueChanged(null);
                }
            } catch (Exception e) {
                // do nothing
            }

            try {
                if (((Datum) currentNode.getUserObject()).getData() instanceof Alignment) {
                    SeqViewerPanel.removeInstance((Alignment) ((Datum) currentNode.getUserObject()).getData());
                }
            } catch (Exception e) {
                // do nothing
            }
        }


    }

    public Map getDataList() {
        return Collections.unmodifiableMap(myDataSetList);
    }

    public JTree getTree() {
        return myTree;
    }

    /**
     * Returns Tassel data set after complete.
     *
     * @param event event
     */
    public void dataSetReturned(PluginEvent event) {
        DataSet tds = (DataSet) event.getSource();
        if (tds != null) {
            addDataSet(tds, DataTreePanel.NODE_TYPE_DEFAULT);
        }
    }

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    public void progress(PluginEvent event) {
        //do nothing
    }

    class MyTreeModelListener implements javax.swing.event.TreeModelListener {

        DataTreePanel theDTP;

        public MyTreeModelListener(DataTreePanel theDTP) {
            this.theDTP = theDTP;
        }

        public void treeNodesChanged(TreeModelEvent e) {
            DefaultMutableTreeNode node;
            node = (DefaultMutableTreeNode) (e.getTreePath().getLastPathComponent());

            /*
             * If the event lists children, then the changed
             * node is the child of the node we've already
             * gotten.  Otherwise, the changed node and the
             * specified node are the same.
             */
            try {
                int index = e.getChildIndices()[0];
                node = (DefaultMutableTreeNode) (node.getChildAt(index));
            } catch (NullPointerException exc) {
            }
            String s = (String) node.getUserObject();
            theDTP.myLastBookSelected.setName(s);
            node.setUserObject(theDTP.myLastBookSelected);
            System.out.println("The user has finished editing the node.");
            System.out.println("New value: " + node.getUserObject());
        }

        public void treeNodesInserted(TreeModelEvent e) {
        }

        public void treeNodesRemoved(TreeModelEvent e) {
        }

        public void treeStructureChanged(TreeModelEvent e) {
        }
    }
}
