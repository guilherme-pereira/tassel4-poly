/*
 * LinkageDiseqDisplayPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.gui.LinkageDisequilibriumComponent;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;
import javax.swing.event.ChangeEvent;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;

import java.net.URL;
import java.util.List;
import net.maizegenetics.pal.gui.LinkageDisequilibriumMinimapComponent;


/**
 *
 * @author Ed Buckler
 */
public class LinkageDiseqDisplayPlugin extends AbstractDisplayPlugin {

    private int upperCorner = LinkageDisequilibriumComponent.P_VALUE;
    private int lowerCorner = LinkageDisequilibriumComponent.RSQUARE;
    private boolean blockSchematic = true;
    private boolean chromosomalView = false;
    private boolean includeLabels = true;

    /** Creates a new instance of LinkageDiseqDisplayPlugin */
    public LinkageDiseqDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> LDInList = input.getDataOfType(LinkageDisequilibrium.class);
            if (LDInList.size() != 1) {
                String message = "Invalid selection.  Please select one LD result.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }
            LinkageDisequilibrium theLD = (LinkageDisequilibrium) LDInList.get(0).getData();
            if (isInteractive()) {
                try {
                    LinkageDiseqDisplayDialog myDialog = new LinkageDiseqDisplayDialog(this, theLD);
                    myDialog.setLocationRelativeTo(getParentFrame());
                    myDialog.setVisible(true);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create LD plot " + ex);
                } catch (Error er) {
                    er.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create LD plot " + er);
                }
            } else if (getSaveFile() != null) {
//                LinkageDisequilibriumComponent ldc = new LinkageDisequilibriumComponent(theLD, blockSchematic, chromosomalView, theLD.getSiteCount(), Math.max(1, theLD.getSiteCount()), Math.max(1, theLD.getSiteCount()));
                LinkageDisequilibriumComponent ldc = new LinkageDisequilibriumComponent(theLD, blockSchematic, chromosomalView, theLD.getSiteCount(), (int)Math.ceil(theLD.getSiteCount()/2.0), (int)Math.ceil(theLD.getSiteCount()/2.0));
                ldc.setUpperCorner(upperCorner);
                ldc.setLowerCorner(lowerCorner);
                ldc.setSize(getImageWidth(), getImageHeight());
                ldc.setShowLabels(includeLabels);
                saveDataToFile(ldc, getSaveFile());
            }

            return null;
        } finally {
            fireProgress(100);
        }
    }

    public int getUpperCorner() {
        return upperCorner;
    }

    public void setUpperCorner(int upperCorner) {
        this.upperCorner = upperCorner;
    }

    public int getLowerCorner() {
        return lowerCorner;
    }

    public void setLowerCorner(int lowerCorner) {
        this.lowerCorner = lowerCorner;
    }

    public boolean isBlockSchematic() {
        return blockSchematic;
    }

    public void setBlockSchematic(boolean blockSchematic) {
        this.blockSchematic = blockSchematic;
    }

    public void setShowLabels(boolean includeLabels) {
        this.includeLabels = includeLabels;
    }

    public boolean isChromosomalView() {
        return chromosomalView;
    }

    public void setChromosomalView(boolean chromosomalView) {
        this.chromosomalView = chromosomalView;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = LinkageDiseqDisplayPlugin.class.getResource("images/LDPlot.gif");
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
        return "LD Plot";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Display Linkage Disequilibrium";
    }
}

class LinkageDiseqDisplayDialog extends JDialog {

    final int myDefaultViewableSize = 35;

    int maxWindowSize;
    int myMaxWindowSize = 1000;

    JPanel panel1 = new JPanel();
    JButton okButton = new JButton();
    JPanel linkPanel = new JPanel();
    JPanel mapPanel = new JPanel();
    LinkageDisequilibriumComponent ldFigurePanel;
    LinkageDisequilibriumMinimapComponent ldMapPanel;
    BorderLayout borderLayout1 = new BorderLayout();
    LinkageDisequilibrium theLinkageDisequilibrium;
    LinkageDiseqDisplayPlugin theLinkageDiseqDisplayPlugin;
    BorderLayout borderLayout2 = new BorderLayout();
    JPanel jPanel1 = new JPanel();
    JButton saveButton = new JButton();
    JButton saveAllButton = new JButton();

    JLabel upperSqrLabel;
    JLabel lowerSqrLabel;
    JComboBox upperSqrSelector;
    JComboBox lowerSqrSelector;

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JCheckBox schematicCheckBox = new JCheckBox();

    //locks scroll bars together
    JCheckBox myLockScrollBarsCheckBox = new JCheckBox();
    boolean myLockedBars = false;

    JScrollBar verticalScrollBar = new JScrollBar(Adjustable.VERTICAL);
    JScrollBar horizontalScrollBar = new JScrollBar(Adjustable.HORIZONTAL);
//  JComboBox formatComboBox = new JComboBox();

//  Windowed Viewer
    JPanel viewerLocationOptionsPanel = new JPanel();

  //window size
    JLabel windowSizeLabel = new JLabel();
    JSlider windowSizeSlider;
    JTextField windowSizeText = new JTextField();

  //x axis
    JSlider windowXSlider;


//  Y axis panel
    JPanel yAxisPanel = new JPanel();

  //y axis
    JSlider windowYSlider;

    // remembers current location
    int myXPos;
    int myYPos;

    public LinkageDiseqDisplayDialog(LinkageDiseqDisplayPlugin theQAF, LinkageDisequilibrium theLinkageDisequilibrium) {
        super(theQAF.getParentFrame(), "Linkage Disequilibrium", false);

        this.theLinkageDiseqDisplayPlugin = theQAF;
        this.theLinkageDisequilibrium = theLinkageDisequilibrium;
        maxWindowSize = theLinkageDisequilibrium.getSiteCount() + theLinkageDisequilibrium.getAlignment().getNumLoci() - 1;

        myXPos = (int)Math.ceil(maxWindowSize/2.0);
        myYPos = (int)Math.ceil(maxWindowSize/2.0);

        try {
            jbInit();
            if (theLinkageDisequilibrium.getSiteCount() > myDefaultViewableSize) {
                ldFigurePanel = new LinkageDisequilibriumComponent(theLinkageDisequilibrium, true, false, myDefaultViewableSize, theLinkageDisequilibrium.getSiteCount()/2, theLinkageDisequilibrium.getSiteCount()/2);
            } else {
                ldFigurePanel = new LinkageDisequilibriumComponent(theLinkageDisequilibrium, true, false, maxWindowSize, (int)Math.ceil(maxWindowSize/2.0), (int)Math.ceil(maxWindowSize/2.0));
            }
            linkPanel.add(ldFigurePanel, BorderLayout.CENTER);
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
            JOptionPane.showMessageDialog(this.getParent(), "Unable to create LD plot " + ex);
        }
        repaint();
    }

    void jbInit() throws Exception {

        jPanel1.setLayout(gridBagLayout1);
        jPanel1.setToolTipText("");

//  Window Viewer init

        if (theLinkageDisequilibrium.getSiteCount() > myDefaultViewableSize) {

            viewerLocationOptionsPanel.setLayout(new GridBagLayout());
            yAxisPanel.setLayout(new GridBagLayout());

            //Window size
            windowSizeLabel.setText("Select viewable size");
            windowSizeSlider = new JSlider(1, Math.min(maxWindowSize, myMaxWindowSize), myDefaultViewableSize);
            windowSizeText.setText(String.valueOf(windowSizeSlider.getValue()));

            windowSizeSlider.addMouseListener(new java.awt.event.MouseListener() {

                public void mouseClicked(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }

                public void mousePressed(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }

                public void mouseReleased(MouseEvent me) {
                    sizeSlider_actionPerformed(me);
                }

                public void mouseEntered(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }

                public void mouseExited(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }
            });

            windowSizeText.addKeyListener(new java.awt.event.KeyListener() {

                public void keyTyped(KeyEvent ke) {
    //                countTextField_keyTyped(ke);
                }

                public void keyPressed(KeyEvent ke) {
    //                throw new UnsupportedOperationException("Not supported yet.");
                }

                public void keyReleased(KeyEvent ke) {
                    sizeTextField_keyTyped(ke);
                }
            });

            //X coords
//            windowXSlider = new JSlider(myDefaultViewableSize/2, numSites-(myDefaultViewableSize/2+(myDefaultViewableSize%2)+numChroms-1), numSites/2);
            windowXSlider = new JSlider((int)Math.ceil(myDefaultViewableSize/2.0), (int)Math.ceil(maxWindowSize-myDefaultViewableSize/2.0), (int)Math.ceil(maxWindowSize/2.0));

            windowXSlider.addChangeListener(new javax.swing.event.ChangeListener() {

                public void stateChanged(ChangeEvent ce) {
                    xSlider_actionPerformed(ce);
                }
//
            });

            //Y coords
//            windowYSlider = new JSlider(myDefaultViewableSize/2, numSites-(myDefaultViewableSize/2+(myDefaultViewableSize%2)+numChroms-1), numSites/2);
            windowYSlider = new JSlider((int)Math.ceil(myDefaultViewableSize/2.0), (int)Math.ceil(maxWindowSize-myDefaultViewableSize/2.0), (int)Math.ceil(maxWindowSize/2.0));
            windowYSlider.setOrientation(JSlider.VERTICAL);
            windowYSlider.setInverted(true);

            windowYSlider.addChangeListener(new javax.swing.event.ChangeListener() {

                public void stateChanged(ChangeEvent ce) {
                    ySlider_actionPerformed(ce);
                }

            });

            //Minimap
            ldMapPanel = new LinkageDisequilibriumMinimapComponent(theLinkageDisequilibrium, myDefaultViewableSize, (int)Math.ceil(maxWindowSize/2.0), (int)Math.ceil(maxWindowSize/2.0));
            mapPanel.add(ldMapPanel, BorderLayout.CENTER);
            mapPanel.setBorder(BorderFactory.createEtchedBorder());

            ldMapPanel.addMouseListener(new java.awt.event.MouseListener() {

                public void mouseClicked(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }

                public void mousePressed(MouseEvent me) {
                    mapPanel_mouseEvent(me);
                }

                public void mouseReleased(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }

                public void mouseEntered(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }

                public void mouseExited(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }
            });

            ldMapPanel.addMouseMotionListener(new java.awt.event.MouseMotionListener() {

                public void mouseDragged(MouseEvent me) {
                    mapPanel_mouseEvent(me);
                }

                public void mouseMoved(MouseEvent me) {
//                    throw new UnsupportedOperationException("Not supported yet.");
                }
            });

            //Add window size to panel
            viewerLocationOptionsPanel.add(windowSizeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.LAST_LINE_END, GridBagConstraints.NONE, new Insets(0, 5, 8, 5), 0, 0));
            viewerLocationOptionsPanel.add(windowSizeSlider, new GridBagConstraints(1, 0, 2, 1, 0.9, 1.0, GridBagConstraints.PAGE_END, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
            viewerLocationOptionsPanel.add(windowSizeText, new GridBagConstraints(3, 0, 1, 1, 0.1, 1.0, GridBagConstraints.LAST_LINE_START, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 5), 0, 0));
            viewerLocationOptionsPanel.add(mapPanel, new GridBagConstraints(4, 0, 1, 1, 0.2, 0.2, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));

            //Add X coords to panel
            jPanel1.add(windowXSlider, new GridBagConstraints(0, 0, 6, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 172, 0, 85), 0, 0));

            //Add Y coords to panel
            yAxisPanel.add(windowYSlider, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(10, 0, 30, 0), 0, 0));

            getContentPane().add(viewerLocationOptionsPanel, BorderLayout.NORTH);
            getContentPane().add(mapPanel, BorderLayout.LINE_START);
            getContentPane().add(yAxisPanel, BorderLayout.LINE_END);
        }

        String[] valueOptions = {"P-Value", "R Squared", "D Prime"};

        upperSqrLabel = new JLabel("Upper triangle:");
        upperSqrSelector = new JComboBox(valueOptions);
        upperSqrSelector.setSelectedIndex(1);
        upperSqrSelector.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                upperSqrSelector_actionPerformed(e);
            }
        });

        lowerSqrLabel = new JLabel("Lower triangle:");
        lowerSqrSelector = new JComboBox(valueOptions);
        lowerSqrSelector.setSelectedIndex(0);
        lowerSqrSelector.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                lowerSqrSelector_actionPerformed(e);
            }
        });

        panel1.setLayout(borderLayout2);
        okButton.setText("Close");
        okButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });
        linkPanel.setBorder(BorderFactory.createEtchedBorder());
        linkPanel.setLayout(borderLayout1);
        panel1.setPreferredSize(new Dimension(600, 600));

        saveButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveButton_actionPerformed(e);
            }
        });
        saveButton.setText("Save");

        saveAllButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveAllButton_actionPerformed(e);
            }
        });
        saveAllButton.setText("Save All");

        schematicCheckBox.setSelected(true);
        schematicCheckBox.setText("Schematic");
        schematicCheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                schematicCheckBox_actionPerformed(e);
            }
        });

        myLockScrollBarsCheckBox.setSelected(false);
        myLockScrollBarsCheckBox.setText("Lock Y Axis to X");
        myLockScrollBarsCheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                lockScrollBarsCheckBox_actionPerformed(e);
            }
        });

        getContentPane().add(panel1, BorderLayout.CENTER);
        panel1.add(linkPanel, BorderLayout.CENTER);
        this.getContentPane().add(jPanel1, BorderLayout.SOUTH);
        jPanel1.add(upperSqrLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 10, 0, 0), 5, 0));
        jPanel1.add(lowerSqrLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 10, 0, 0), 5, 0));
        jPanel1.add(upperSqrSelector, new GridBagConstraints(1, 1, 1, 1, 0.5, 0.5, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(lowerSqrSelector, new GridBagConstraints(1, 2, 1, 1, 0.5, 0.5, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(okButton, new GridBagConstraints(5, 2, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 0, 0, 15), 0, 0));
        jPanel1.add(myLockScrollBarsCheckBox, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(saveButton, new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(saveAllButton, new GridBagConstraints(5, 1, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 0, 0, 15), 0, 0));
        jPanel1.add(schematicCheckBox, new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    }

    void okButton_actionPerformed(ActionEvent e) {
        dispose();
    }/////////end calculateDisequilibrium

    void saveButton_actionPerformed(ActionEvent e) {
        //     String s=formatComboBox.getSelectedItem().toString();
        this.theLinkageDiseqDisplayPlugin.saveDataToFile(ldFigurePanel);
    }

    void saveAllButton_actionPerformed(ActionEvent e) {

        LinkageDisequilibriumComponent newGraph = new LinkageDisequilibriumComponent(theLinkageDisequilibrium, true, false, maxWindowSize, (int)Math.ceil(maxWindowSize/2.0), (int)Math.ceil(maxWindowSize/2.0));

        if (upperSqrSelector.getSelectedIndex() == 0) {
            newGraph.setUpperCorner(LinkageDisequilibriumComponent.P_VALUE);
        } else if (upperSqrSelector.getSelectedIndex() == 1) {
            newGraph.setUpperCorner(LinkageDisequilibriumComponent.RSQUARE);
        } else if (upperSqrSelector.getSelectedIndex() == 2) {
            newGraph.setUpperCorner(LinkageDisequilibriumComponent.DPRIME);
        }

        if (lowerSqrSelector.getSelectedIndex() == 0) {
            newGraph.setLowerCorner(LinkageDisequilibriumComponent.P_VALUE);
        } else if (lowerSqrSelector.getSelectedIndex() == 1) {
            newGraph.setLowerCorner(LinkageDisequilibriumComponent.RSQUARE);
        } else if (lowerSqrSelector.getSelectedIndex() == 2) {
            newGraph.setLowerCorner(LinkageDisequilibriumComponent.DPRIME);
        }

        newGraph.setGraphSize(4096, 4096);
        theLinkageDiseqDisplayPlugin.saveDataToFile(newGraph);
    }

    void upperSqrSelector_actionPerformed(ActionEvent e) {
        if (upperSqrSelector.getSelectedIndex() == 0) {
            ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.P_VALUE);
        } else if (upperSqrSelector.getSelectedIndex() == 1) {
            ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.RSQUARE);
        } else if (upperSqrSelector.getSelectedIndex() == 2) {
            ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.DPRIME);
        }
        repaint();
    }

    void lowerSqrSelector_actionPerformed(ActionEvent e) {
        if (lowerSqrSelector.getSelectedIndex() == 0) {
            ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.P_VALUE);
        } else if (lowerSqrSelector.getSelectedIndex() == 1) {
            ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.RSQUARE);
        } else if (lowerSqrSelector.getSelectedIndex() == 2) {
            ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.DPRIME);
        }
        repaint();
    }

    void upDPrimeRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.DPRIME);
        repaint();
    }

    void upRSqrRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.RSQUARE);
        repaint();
    }

    void upPRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.P_VALUE);
        repaint();
    }

    void lowRSqrRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.RSQUARE);
        repaint();
    }

    void lowPRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.P_VALUE);
        repaint();
    }

    void lowDPrimeRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.DPRIME);
        repaint();
    }

    void geneRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setScaleOfView(false);
        repaint();
    }

    void chromoRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setScaleOfView(true);
        repaint();
    }

    void schematicCheckBox_actionPerformed(ActionEvent e) {
        ldFigurePanel.setShowSchematic(schematicCheckBox.isSelected());
        repaint();
    }

    void lockScrollBarsCheckBox_actionPerformed(ActionEvent e) {
        if ( myLockScrollBarsCheckBox.isSelected() ) {
            myLockedBars = true;
        } else {
            myLockedBars = false;
        }
    }

    void sizeSlider_actionPerformed(MouseEvent me) {
        windowSizeText.setText(String.valueOf(windowSizeSlider.getValue()));

//        windowXSlider.setValue((int)Math.ceil(maxWindowSize/2.0));
        windowXSlider.setMinimum((int)Math.ceil(windowSizeSlider.getValue()/2.0));
        windowXSlider.setMaximum((int)Math.ceil(maxWindowSize-windowSizeSlider.getValue()/2.0));
        windowXSlider.setValue(Math.min(windowXSlider.getMaximum(), Math.max(windowXSlider.getMinimum(), myXPos)));

//        windowYSlider.setValue((int)Math.ceil(maxWindowSize/2.0));
        windowYSlider.setMinimum((int)Math.ceil(windowSizeSlider.getValue()/2.0));
        windowYSlider.setMaximum((int)Math.ceil(maxWindowSize-windowSizeSlider.getValue()/2.0));
        windowYSlider.setValue(Math.min(windowYSlider.getMaximum(), Math.max(windowYSlider.getMinimum(), myYPos)));

        int ldMeasureLower = 0;
        int ldMeasureUpper = 0;

        if (lowerSqrSelector.getSelectedIndex() == 0) {
            ldMeasureLower = LinkageDisequilibriumComponent.P_VALUE;
        } else if (lowerSqrSelector.getSelectedIndex() == 1) {
            ldMeasureLower = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureLower = LinkageDisequilibriumComponent.DPRIME;
        }

        if (upperSqrSelector.getSelectedIndex() == 0) {
            ldMeasureUpper = LinkageDisequilibriumComponent.P_VALUE;
        } else if (upperSqrSelector.getSelectedIndex() == 1) {
            ldMeasureUpper = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureUpper = LinkageDisequilibriumComponent.DPRIME;
        }

        ldFigurePanel.setWindowSize(windowSizeSlider.getValue(), windowXSlider.getValue(), windowYSlider.getValue(), ldMeasureLower, ldMeasureUpper);
        ldMapPanel.setWindowSize(windowSizeSlider.getValue(), windowXSlider.getValue(), windowYSlider.getValue());

        repaint();
    }

    void sizeTextField_keyTyped(KeyEvent ke) {
        try {
            if (!windowSizeText.getText().equals("")) {
                int value = Integer.valueOf(windowSizeText.getText());
                if (value >= windowSizeSlider.getMinimum() && value <= windowSizeSlider.getMaximum()) {
                    windowSizeSlider.setValue(value);
                } else if (value <= windowSizeSlider.getMinimum()) {
                    windowSizeSlider.setValue(windowSizeSlider.getMinimum());
                    windowSizeText.setText(String.valueOf(windowSizeSlider.getMinimum()));
                } else if (value >= windowSizeSlider.getMaximum()) {
                    windowSizeSlider.setValue(windowSizeSlider.getMaximum());
                    windowSizeText.setText(String.valueOf(windowSizeSlider.getMaximum()));
                }
            }
        } catch (NumberFormatException nfe) {
            windowSizeText.setText(String.valueOf(windowSizeSlider.getValue()));
        }
    }

    void xSlider_actionPerformed(ChangeEvent ce) {
        int ldMeasureLower = 0;
        int ldMeasureUpper = 0;

        if (lowerSqrSelector.getSelectedIndex() == 0) {
            ldMeasureLower = LinkageDisequilibriumComponent.P_VALUE;
        } else if (lowerSqrSelector.getSelectedIndex() == 1) {
            ldMeasureLower = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureLower = LinkageDisequilibriumComponent.DPRIME;
        }

        if (upperSqrSelector.getSelectedIndex() == 0) {
            ldMeasureUpper = LinkageDisequilibriumComponent.P_VALUE;
        } else if (upperSqrSelector.getSelectedIndex() == 1) {
            ldMeasureUpper = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureUpper = LinkageDisequilibriumComponent.DPRIME;
        }

        ldFigurePanel.setWindowX(windowXSlider.getValue(), ldMeasureLower, ldMeasureUpper);
        ldMapPanel.setWindowX(windowXSlider.getValue());

        if ( myLockedBars ) {
            int diff = myXPos - windowXSlider.getValue();
            System.out.println(diff + "");
            if ( diff < 0) {
                windowYSlider.setValue(Math.max(windowYSlider.getMinimum(), windowYSlider.getValue() - diff));
            } else {
                windowYSlider.setValue(Math.min(windowYSlider.getMaximum(), windowYSlider.getValue() - diff));
            }

            ldFigurePanel.setWindowY(windowYSlider.getValue(), ldMeasureLower, ldMeasureUpper);
            ldMapPanel.setWindowY(windowYSlider.getValue());
            myYPos = windowYSlider.getValue();
        }

        myXPos = windowXSlider.getValue();

        repaint();
    }

    void ySlider_actionPerformed(ChangeEvent ce) {
        int ldMeasureLower = 0;
        int ldMeasureUpper = 0;

        if (lowerSqrSelector.getSelectedIndex() == 0) {
            ldMeasureLower = LinkageDisequilibriumComponent.P_VALUE;
        } else if (lowerSqrSelector.getSelectedIndex() == 1) {
            ldMeasureLower = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureLower = LinkageDisequilibriumComponent.DPRIME;
        }

        if (upperSqrSelector.getSelectedIndex() == 0) {
            ldMeasureUpper = LinkageDisequilibriumComponent.P_VALUE;
        } else if (upperSqrSelector.getSelectedIndex() == 1) {
            ldMeasureUpper = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureUpper = LinkageDisequilibriumComponent.DPRIME;
        }

        ldFigurePanel.setWindowY(windowYSlider.getValue(), ldMeasureLower, ldMeasureUpper);
        ldMapPanel.setWindowY(windowYSlider.getValue());

        myYPos = windowYSlider.getValue();

        repaint();
    }

    void mapPanel_mouseEvent(MouseEvent me) {

        double mouseX = me.getLocationOnScreen().getX();
        double mouseY = me.getLocationOnScreen().getY();

        double panelX = ldMapPanel.getLocationOnScreen().getX();
        double panelY = ldMapPanel.getLocationOnScreen().getY();

        int mapSize = ldMapPanel.getMapSize();

        double newX = (mouseX-panelX)/mapSize;
        double newY = (mouseY-panelY)/mapSize;

        int newXPos = Math.min(Math.max((int)(maxWindowSize*newX), windowXSlider.getMinimum()), windowXSlider.getMaximum());
        int newYPos = Math.min(Math.max((int)(maxWindowSize*newY), windowYSlider.getMinimum()), windowYSlider.getMaximum());

        int ldMeasureLower = 0;
        int ldMeasureUpper = 0;

        if (lowerSqrSelector.getSelectedIndex() == 0) {
            ldMeasureLower = LinkageDisequilibriumComponent.P_VALUE;
        } else if (lowerSqrSelector.getSelectedIndex() == 1) {
            ldMeasureLower = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureLower = LinkageDisequilibriumComponent.DPRIME;
        }

        if (upperSqrSelector.getSelectedIndex() == 0) {
            ldMeasureUpper = LinkageDisequilibriumComponent.P_VALUE;
        } else if (upperSqrSelector.getSelectedIndex() == 1) {
            ldMeasureUpper = LinkageDisequilibriumComponent.RSQUARE;
        } else {
            ldMeasureUpper = LinkageDisequilibriumComponent.DPRIME;
        }

        windowXSlider.setValue(newXPos);
        windowYSlider.setValue(newYPos);

        ldFigurePanel.setWindowX(newXPos, ldMeasureLower, ldMeasureUpper);
        ldFigurePanel.setWindowY(newYPos, ldMeasureLower, ldMeasureUpper);

        ldMapPanel.setWindowX(newXPos);
        ldMapPanel.setWindowY(newYPos);

        myXPos = windowXSlider.getValue();
        myYPos = windowYSlider.getValue();

        repaint();
    }

    public int getWindowSizeSelection() {
        return windowSizeSlider.getValue();
    }

    public int getWindowXSelection() {
        return windowXSlider.getValue();
    }

    public int getWindowYSelection() {
        return windowYSlider.getValue();
    }
}