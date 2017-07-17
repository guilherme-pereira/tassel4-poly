package net.maizegenetics.tassel;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;

import net.maizegenetics.prefs.TasselPrefs;

/**
 * @author terryc
 */
public class PreferencesDialog extends JDialog {
    
    private final static int TEXT_FIELD_WIDTH = 10;
    private final static Font HEADING_FONT = new Font(null, Font.BOLD, 14);
    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JRadioButton myStrictButton = new JRadioButton("Strict");
    private JRadioButton myNonStrictButton = new JRadioButton("Non-Strict");
    private JRadioButton myFirstLevelButton = new JRadioButton("First Level");
    private JTextField myMaxRetainAlleles = new JTextField(TEXT_FIELD_WIDTH);
    private JCheckBox myRetainRareAlleles = new JCheckBox("Retain Rare Alleles");
    
    public PreferencesDialog() {
        super((Frame) null, "Preferences...", true);
        try {
            init();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private void init() throws Exception {
        
        setTitle("Preferences...");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);
        
        Container contentPane = getContentPane();
        
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.LEFT_ALIGNMENT);
        
        result.add(getStrictPanel());
        
        result.add(Box.createRigidArea(new Dimension(1, 30)));
        
        result.add(getGenotypeStoringPanel());
        
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        
        result.add(getButtons());
        
        result.add(Box.createRigidArea(new Dimension(1, 20)));
        
        contentPane.add(result);
        
        pack();
        
        setResizable(false);
        
    }
    
    private JPanel getGenotypeStoringPanel() {
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        
        JLabel alignPrefsLabel = new JLabel("Alignment Preferences...");
        alignPrefsLabel.setFont(HEADING_FONT);
        result.add(alignPrefsLabel);
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        
        myMaxRetainAlleles.setText(String.valueOf(TasselPrefs.getAlignmentMaxAllelesToRetain()));
        result.add(getLine("Max Alleles Retained: ", myMaxRetainAlleles));
        
        myRetainRareAlleles.setSelected(TasselPrefs.getAlignmentRetainRareAlleles());
        result.add(myRetainRareAlleles);
        
        return result;
    }
    
    private JPanel getLine(String label, JTextField ref) {
        
        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));
        
        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);
        
        return result;
        
    }
    
    private JPanel getStrictPanel() {
        
        JPanel inputs = new JPanel();
        BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
        inputs.setLayout(layout);
        inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        
        inputs.add(Box.createRigidArea(new Dimension(1, 10)));
        
        inputs.add(getStrictLabel());
        
        inputs.add(getStrictOptionPanel());
        
        return inputs;
        
    }
    
    private JPanel getStrictLabel() {
        
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        
        int borderWidth = 15;
        
        String msg1 = "Choose Strict or Non-Strict...  If Strict, taxa names";
        String msg2 = "must match exactly. If Non-Strict, taxa names must";
        String msg3 = "match up to least descriptive taxa.";
        JLabel jmsg1 = new JLabel(msg1, JLabel.CENTER);
        jmsg1.setFont(HEADING_FONT);
        jmsg1.setBorder(BorderFactory.createEmptyBorder(0, borderWidth, 0, borderWidth));
        JLabel jmsg2 = new JLabel(msg2, JLabel.CENTER);
        jmsg2.setFont(HEADING_FONT);
        jmsg2.setBorder(BorderFactory.createEmptyBorder(0, borderWidth, 0, borderWidth));
        JLabel jmsg3 = new JLabel(msg3, JLabel.CENTER);
        jmsg3.setFont(HEADING_FONT);
        jmsg3.setBorder(BorderFactory.createEmptyBorder(0, borderWidth, 0, borderWidth));
        
        result.add(Box.createRigidArea(new Dimension(1, 20)));
        result.add(jmsg1);
        result.add(jmsg2);
        result.add(jmsg3);
        result.add(Box.createRigidArea(new Dimension(1, 20)));

        //jLabel1.setFont(new Font("Dialog", Font.BOLD, 16));
        return result;
        
    }
    
    private JPanel getStrictOptionPanel() {
        
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        result.setBorder(BorderFactory.createEtchedBorder());
        
        myButtonGroup.add(myStrictButton);
        myButtonGroup.add(myNonStrictButton);
        myButtonGroup.add(myFirstLevelButton);
        
        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();
        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
            myStrictButton.setSelected(true);
        } else if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NonStrict) {
            myNonStrictButton.setSelected(true);
        } else if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NumLevels) {
            myFirstLevelButton.setSelected(true);
        } else {
            throw new IllegalStateException("PreferencesDialog: getStrictOptionPanel: Unknown Taxa Join Type: " + type.toString());
        }
        
        result.add(myNonStrictButton);
        result.add(myStrictButton);
        result.add(myFirstLevelButton);

        // result.add(Box.createRigidArea(new Dimension(1, 20)));

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
    
    private void okButton_actionPerformed(ActionEvent e) {
        
        if (myStrictButton.isSelected()) {
            TasselPrefs.putIDJoinStrict(TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict);
        } else if (myFirstLevelButton.isSelected()) {
            TasselPrefs.putIDJoinStrict(TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NumLevels);
        } else {
            TasselPrefs.putIDJoinStrict(TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NonStrict);
        }
        
        try {
            int maxNumAlleles = Integer.valueOf(myMaxRetainAlleles.getText());
            if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
                JOptionPane.showMessageDialog(this, "Max Alleles Retained must be between 1 and 14 inclusive", "Error", JOptionPane.ERROR_MESSAGE);
                return;
            }
            TasselPrefs.putAlignmentMaxAllelesToRetain(maxNumAlleles);
        } catch (NumberFormatException ex1) {
            JOptionPane.showMessageDialog(this, "Max Alleles Retained must be integer number", "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        
        TasselPrefs.putAlignmentRetainRareAlleles(myRetainRareAlleles.isSelected());
        
        setVisible(false);
    }
    
    private void cancelButton_actionPerformed(ActionEvent e) {
        
        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();
        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
            myStrictButton.setSelected(true);
        } else if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NumLevels) {
            myFirstLevelButton.setSelected(true);
        } else {
            myNonStrictButton.setSelected(true);
        }
        
        myMaxRetainAlleles.setText(String.valueOf(TasselPrefs.getAlignmentMaxAllelesToRetain()));
        myRetainRareAlleles.setSelected(TasselPrefs.getAlignmentRetainRareAlleles());
        
        setVisible(false);
    }
}
