package net.maizegenetics.wizard.panels;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

import net.maizegenetics.wizard.*;
import java.io.File;
import javax.swing.event.ChangeListener;

import net.maizegenetics.baseplugins.FileLoadPlugin.TasselFileType;


public class TestPanel2 extends JPanel {
 
    private javax.swing.JLabel anotherBlankSpace;
    private javax.swing.JLabel blankSpace;
    private javax.swing.JLabel blankSpace2;
    private javax.swing.JLabel blankSpace3;
    private javax.swing.ButtonGroup connectorGroup;
//    private javax.swing.JRadioButton ethernetRJRadioButton;
//    private javax.swing.JRadioButton ethernetTenRadioButton;
//    private javax.swing.JCheckBox jCheckBox1;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel browsePanel;
//    private javax.swing.JRadioButton notInventedYetRadioButton;
//    private javax.swing.JRadioButton serialParallelRadioButton;
    private javax.swing.JLabel welcomeTitle;
//    private javax.swing.JRadioButton wirelessRadioButton;
    private javax.swing.JLabel yetAnotherBlankSpace1;
    private javax.swing.JRadioButton hapmapButton = new JRadioButton("Load Hapmap");
    private javax.swing.JRadioButton plinkButton = new JRadioButton("Load Plink");
    private javax.swing.JRadioButton seqAlignButton = new JRadioButton("Load Sequence Alignment");
    private javax.swing.JRadioButton fastaButton = new JRadioButton("Load FASTA file");
    private javax.swing.JRadioButton polyAlignButton = new JRadioButton("Load Polymorphism Alignment (custom)");
    private javax.swing.JRadioButton annoAlignButton = new JRadioButton("Load Annotated Alignment (custom)");
    private javax.swing.JRadioButton numCovButton = new JRadioButton("Load Numerial Trait Data or Covariates");
    private javax.swing.JRadioButton sqrMatButton = new JRadioButton("Load Square Numerical Matrix (ex. kinship) (phylip)");
    private javax.swing.JRadioButton guessButton = new JRadioButton("I will make my best guess");
    private javax.swing.JButton browseButton = new JButton("Browse...");
    private javax.swing.JTextField browseTextField = new JTextField();

    private javax.swing.JFileChooser fileChooser = new JFileChooser();
    
    private JPanel contentPanel;
    private JLabel iconLabel;
    private JSeparator separator;
    private JLabel textLabel;
    private JPanel titlePanel;

    private JLabel jLabel2;
    private JLabel jLabel3;
    private JLabel jLabel4;
    private JLabel jLabel5;
    private JLabel jLabel6;
    private JLabel jLabel7;
    private JLabel jLabel8;
    private JLabel jLabel9;

    public TasselFileType fileType = TasselFileType.Unknown;
        
    public TestPanel2() {
     
        super();
                
        contentPanel = getContentPanel();
        contentPanel.setBorder(new EmptyBorder(new Insets(10, 10, 10, 10)));

        ImageIcon icon = getImageIcon();
        
        titlePanel = new javax.swing.JPanel();
        textLabel = new javax.swing.JLabel();
        iconLabel = new javax.swing.JLabel();
        separator = new javax.swing.JSeparator();

        setLayout(new java.awt.BorderLayout());


        titlePanel.setLayout(new java.awt.BorderLayout());
        //titlePanel.setBackground(Color.gray);
        
        //textLabel.setBackground(Color.gray);
        textLabel.setFont(new Font("MS Sans Serif", Font.BOLD, 14));
        textLabel.setText("Select Data File");
        textLabel.setBorder(new EmptyBorder(new Insets(10, 10, 10, 10)));
        textLabel.setOpaque(true);

        iconLabel.setBackground(Color.gray);
        if (icon != null)
            iconLabel.setIcon(icon);
        
        titlePanel.add(textLabel, BorderLayout.CENTER);
        titlePanel.add(iconLabel, BorderLayout.EAST);
        titlePanel.add(separator, BorderLayout.SOUTH);
        
        add(titlePanel, BorderLayout.NORTH);
        JPanel secondaryPanel = new JPanel();
        secondaryPanel.add(contentPanel, BorderLayout.NORTH);
        add(secondaryPanel, BorderLayout.WEST);

    }  
    
//    public void addCheckBoxActionListener(ActionListener l) {
//        jCheckBox1.addActionListener(l);
//    }
    
//    public boolean isCheckBoxSelected() {
//        return jCheckBox1.isSelected();
//    }

//    public void addTextFieldActionListener(ActionListener l) {
//        browseTextField.addActionListener(l);
//    }
//
//    public void addBrowseButtonActionListener(ActionListener l) {
//        browseButton.addActionListener(l);
//    }

//    public void addFileChooserMouseListener(MouseListener l) {
//        fileChooser.addMouseListener(l);
//    }

    public void addPanelMouseListener(MouseListener l) {
        contentPanel.addMouseListener(l);
    }

    public void addTextFieldKeyListener(KeyListener l) {
        browseTextField.addKeyListener(l);
    }

    public void addTextFieldMouseListener(MouseListener l) {
        browseTextField.addMouseListener(l);
    }

    public void addButtonGroupMouseListener(MouseListener l) {
        hapmapButton.addMouseListener(l);
        plinkButton.addMouseListener(l);
        seqAlignButton.addMouseListener(l);
        fastaButton.addMouseListener(l);
        polyAlignButton.addMouseListener(l);
        annoAlignButton.addMouseListener(l);
        numCovButton.addMouseListener(l);
        sqrMatButton.addMouseListener(l);
        guessButton.addMouseListener(l);
    }

    public boolean isTextFieldFilled() {
//        System.out.println(browseTextField.getText().length() == 0);
        return !(browseTextField.getText().length() == 0);
    }
    
    public String getRadioButtonSelected() {
        return connectorGroup.getSelection().getActionCommand();
    }

    public String getBrowseFieldText() {
        return browseTextField.getText();
    }

    public TasselFileType getFileType() {
        return fileType;
    }
    
    private JPanel getContentPanel() {     
        
        JPanel contentPanel1 = new JPanel();
        
        connectorGroup = new javax.swing.ButtonGroup();
        welcomeTitle = new javax.swing.JLabel();
        jPanel1 = new javax.swing.JPanel();
        browsePanel = new javax.swing.JPanel();
        blankSpace = new javax.swing.JLabel();
//        wirelessRadioButton = new javax.swing.JRadioButton();
//        ethernetRJRadioButton = new javax.swing.JRadioButton();
//        ethernetTenRadioButton = new javax.swing.JRadioButton();
//        serialParallelRadioButton = new javax.swing.JRadioButton();
//        notInventedYetRadioButton = new javax.swing.JRadioButton();
        
//        wirelessRadioButton.setActionCommand("Wireless Radio");
//        ethernetRJRadioButton.setActionCommand("Ethernet RJ-45");
//        ethernetTenRadioButton.setActionCommand("Ethernet 10 base T");
//        serialParallelRadioButton.setActionCommand("Serial/Parallel");
//        notInventedYetRadioButton.setActionCommand("Not Yet Invented");

        hapmapButton.setActionCommand("Hapmap");
        hapmapButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                hapmapButton_actionPerformed(e);
            }
        });
        plinkButton.setActionCommand("Plink");
        plinkButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                plinkButton_actionPerformed(e);
            }
        });
        seqAlignButton.setActionCommand("Sequence Alignment");
        seqAlignButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                seqAlignButton_actionPerformed(e);
            }
        });
        fastaButton.setActionCommand("FASTA file");
        fastaButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                fastaButton_actionPerformed(e);
            }
        });
        polyAlignButton.setActionCommand("Polymorphism Alignment");
        polyAlignButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                polyAlignButton_actionPerformed(e);
            }
        });
        annoAlignButton.setActionCommand("Annotated Alignment");
        annoAlignButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                annoAlignButton_actionPerformed(e);
            }
        });
        numCovButton.setActionCommand("Numerical Trait Data or Covariates");
        numCovButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                numCovButton_actionPerformed(e);
            }
        });
        sqrMatButton.setActionCommand("Square Numerical Matrix");
        sqrMatButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                sqrMatButton_actionPerformed(e);
            }
        });
        guessButton.setActionCommand("Guess");
        guessButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                guessButton_actionPerformed(e);
            }
        });

        connectorGroup.add(hapmapButton);
        connectorGroup.add(plinkButton);
        connectorGroup.add(seqAlignButton);
        connectorGroup.add(fastaButton);
        connectorGroup.add(polyAlignButton);
        connectorGroup.add(annoAlignButton);
        connectorGroup.add(numCovButton);
        connectorGroup.add(sqrMatButton);
        connectorGroup.add(guessButton);
        guessButton.setSelected(true);
        
        anotherBlankSpace = new javax.swing.JLabel();
//        jCheckBox1 = new javax.swing.JCheckBox();
        yetAnotherBlankSpace1 = new javax.swing.JLabel();
        jLabel1 = new javax.swing.JLabel();

//        ethernetRJRadioButton.setSelected(true);
        
        contentPanel1.setLayout(new java.awt.BorderLayout());

        welcomeTitle.setText("Select the type of data that you would like to Import:");
        contentPanel1.add(welcomeTitle, java.awt.BorderLayout.NORTH);

        jPanel1.setLayout(new java.awt.GridLayout(0, 1));
        //jPanel1.setLayout(new GridBagLayout());

        //jPanel1.add(blankSpace);

        jPanel1.add(hapmapButton);
        jPanel1.add(plinkButton);
        jPanel1.add(seqAlignButton);
        jPanel1.add(fastaButton);
        jPanel1.add(polyAlignButton);
        jPanel1.add(annoAlignButton);
        jPanel1.add(numCovButton);
        jPanel1.add(sqrMatButton);
        jPanel1.add(guessButton);

//        wirelessRadioButton.setText("802.11 b/g");
//        connectorGroup.add(wirelessRadioButton);
//        jPanel1.add(wirelessRadioButton);
//
//        ethernetRJRadioButton.setText("Ethernet RJ-45");
//        connectorGroup.add(ethernetRJRadioButton);
//        jPanel1.add(ethernetRJRadioButton);
//
//        ethernetTenRadioButton.setText("Ethernet 10 base T");
//        connectorGroup.add(ethernetTenRadioButton);
//        jPanel1.add(ethernetTenRadioButton);
//
//        serialParallelRadioButton.setText("Serial/Parallel");
//        connectorGroup.add(serialParallelRadioButton);
//        jPanel1.add(serialParallelRadioButton);
//
//        notInventedYetRadioButton.setText("Something Not Yet Invented But You're Sure You'll Want It");
//        connectorGroup.add(notInventedYetRadioButton);
//        jPanel1.add(notInventedYetRadioButton);

        jPanel1.add(anotherBlankSpace);

        jLabel2 = new JLabel();
        jLabel2.setText("Browse to data file:");
        jPanel1.add(jLabel2);

        //browsePanel.setLayout(new GridLayout(1, 3));
        browsePanel.setLayout(new GridBagLayout());
        browsePanel.add(browseTextField, new GridBagConstraints(0, 0, 2, 1, 1.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 5), 0, 0));
        browsePanel.add(browseButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 0, 0));

        jPanel1.add(browsePanel);

        browseButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                browseButton_actionPerformed(e);
            }
        });


//        jCheckBox1.setText("I agree to laugh at people who chose options other than mine");
//        jPanel1.add(jCheckBox1);

        jPanel1.add(yetAnotherBlankSpace1);

        contentPanel1.add(jPanel1, java.awt.BorderLayout.CENTER);

        //jLabel1.setText("Note that the 'Next' button is disabled until you check the box above.");
        contentPanel1.add(jLabel1, java.awt.BorderLayout.SOUTH);
        
        return contentPanel1;
    }

    private void browseButton_actionPerformed(ActionEvent e) {
        int returnVal = fileChooser.showOpenDialog(jPanel1);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            browseTextField.setText(file.getPath());
        }
    }
    
    private ImageIcon getImageIcon() {
        
        //  Icon to be placed in the upper right corner.
        
        return null;
    }
    
//    connectorGroup.add(hapmapButton);
//    connectorGroup.add(plinkButton);
//    connectorGroup.add(seqAlignButton);
//    connectorGroup.add(fastaButton);
//    connectorGroup.add(polyAlignButton);
//    connectorGroup.add(annoAlignButton);
//    connectorGroup.add(numCovButton);
//    connectorGroup.add(sqrMatButton);
//    connectorGroup.add(guessButton);

    public void hapmapButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Hapmap;
    }

    public void plinkButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Plink;
    }

    public void seqAlignButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Sequence;
    }

    public void fastaButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Fasta;
    }

    public void polyAlignButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Polymorphism;
    }

    public void annoAlignButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Annotated;
    }

    public void numCovButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Numerical;
    }

    public void sqrMatButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.SqrMatrix;
    }

    public void guessButton_actionPerformed(ActionEvent e) {
        fileType = TasselFileType.Unknown;
    }
}
