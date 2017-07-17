package net.maizegenetics.wizard.panels;

import net.maizegenetics.wizard.WizardPanelDescriptor;
//import net.maizegenetics.wizard.*;

import net.maizegenetics.baseplugins.FileLoadPlugin;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;


public class TestPanel2Descriptor extends WizardPanelDescriptor implements KeyListener, MouseListener {
    
    public static final String IDENTIFIER = "LOAD_FILE_PANEL";
    
    TestPanel2 panel2;
    
    public TestPanel2Descriptor() {
        
        panel2 = new TestPanel2();
//        panel2.addCheckBoxActionListener(this);
        panel2.addTextFieldKeyListener(this);
        panel2.addTextFieldMouseListener(this);
        panel2.addPanelMouseListener(this);
        panel2.addButtonGroupMouseListener(this);
        
        setPanelDescriptorIdentifier(IDENTIFIER);
        setPanelComponent(panel2);
   
    }


    public void aboutToHidePanel() {
        FileLoadPlugin fileLoader = new FileLoadPlugin(null, false);
        fileLoader.setTheFileType(panel2.getFileType());
        String[] fileName = new String[1];
        fileName[0] = panel2.getBrowseFieldText();
        fileLoader.setOpenFiles(fileName);
        data = fileLoader.performFunction(null);
    }

    public void addWizardMouseListener(MouseListener l) {
        this.getWizard().getDialog().addMouseListener(l);
    }
    
    public Object getNextPanelDescriptor() {
        return TestPanel3Descriptor.IDENTIFIER;
    }
    
    public Object getBackPanelDescriptor() {
        return TestPanel1Descriptor.IDENTIFIER;
    }

    public void aboutToDisplayPanel() {
        setNextButtonAccordingToTextField();
    }

    public void keyReleased(KeyEvent e) {
        setNextButtonAccordingToTextField();
    }

    public void keyTyped(KeyEvent e) {
        setNextButtonAccordingToTextField();
    }

    public void keyPressed(KeyEvent e) {
        setNextButtonAccordingToTextField();
    }

//    public void actionPerformed(ActionEvent e) {
//        setNextButtonAccordingToTextField();
//    }

    public void mouseClicked(MouseEvent e) {
        setNextButtonAccordingToTextField();
    }

    public void mousePressed(MouseEvent e) {
        setNextButtonAccordingToTextField();
    }

    public void mouseReleased(MouseEvent e) {
        setNextButtonAccordingToTextField();
    }

    public void mouseEntered(MouseEvent e) {
        setNextButtonAccordingToTextField();
    }

    public void mouseExited(MouseEvent e) {
        setNextButtonAccordingToTextField();
    }

    private void setNextButtonAccordingToTextField() {
         if (panel2.isTextFieldFilled())
            getWizard().setNextFinishButtonEnabled(true);
         else
            getWizard().setNextFinishButtonEnabled(false);

    }
}
