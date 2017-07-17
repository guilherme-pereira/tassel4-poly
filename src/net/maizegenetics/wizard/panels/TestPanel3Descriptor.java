package net.maizegenetics.wizard.panels;

import net.maizegenetics.wizard.WizardPanelDescriptor;
//import net.maizegenetics.wizard.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import net.maizegenetics.baseplugins.LinkageDisequilibriumPlugin;
import net.maizegenetics.baseplugins.GenotypeImputationPlugin;
//import net.maizegenetics.baseplugins.LinkageDiseqDisplayPlugin;
import net.maizegenetics.plugindef.DataSet;


public class TestPanel3Descriptor extends WizardPanelDescriptor implements ActionListener {
    
    public static final String IDENTIFIER = "ANALYSIS_SELECTION_PANEL";
    
    TestPanel3 panel3;
    
    public TestPanel3Descriptor() {
        
        panel3 = new TestPanel3();
        panel3.addCheckBoxActionListener(this);
//        panel3.addLDCheckBoxActionListener(this);
//        panel3.addGLMCheckBoxActionListener(this);
//        panel3.addMLMCheckBoxActionListener(this);
//        panel3.addImputeCheckBoxActionListener(this);
        setPanelDescriptorIdentifier(IDENTIFIER);
        setPanelComponent(panel3);
        
    }

    public Object getNextPanelDescriptor() {
        return TestPanel4Descriptor.IDENTIFIER;
    }
    
    public Object getBackPanelDescriptor() {
        return TestPanel2Descriptor.IDENTIFIER;
    }
    
    
    public void aboutToDisplayPanel() {
        
//        panel3.setProgressValue(0);
//        panel3.setProgressText("Connecting to Server...");

        getWizard().setNextFinishButtonEnabled(false);
        getWizard().setBackButtonEnabled(false);
        
    }
    
//    public void displayingPanel() {
//
//            Thread t = new Thread() {
//
//            public void run() {
//
//                try {
//                    Thread.sleep(2000);
//                    panel3.setProgressValue(25);
//                    panel3.setProgressText("Server Connection Established");
//                    Thread.sleep(500);
//                    panel3.setProgressValue(50);
//                    panel3.setProgressText("Transmitting Data...");
//                    Thread.sleep(3000);
//                    panel3.setProgressValue(75);
//                    panel3.setProgressText("Receiving Acknowledgement...");
//                    Thread.sleep(1000);
//                    panel3.setProgressValue(100);
//                    panel3.setProgressText("Data Successfully Transmitted");
//
//                    getWizard().setNextFinishButtonEnabled(true);
//                    getWizard().setBackButtonEnabled(true);
//
//                } catch (InterruptedException e) {
//
//                    panel3.setProgressValue(0);
//                    panel3.setProgressText("An Error Has Occurred");
//
//                    getWizard().setBackButtonEnabled(true);
//                }
//
//            }
//        };
//
//        t.start();
//    }
    
    public void aboutToHidePanel() {

        if (panel3.isLDChecked()) {
            LinkageDisequilibriumPlugin ld = new LinkageDisequilibriumPlugin(null, false);
            DataSet data2 = ld.performFunction(data);
            this.setData(data2);
        }
        else if (panel3.isImputeChecked()) {
            GenotypeImputationPlugin impute = new GenotypeImputationPlugin(null, false);
            DataSet data2 = impute.performFunction(data);
            this.setData(data2);
        }

//        if (panel3.isLDPlotChecked()) {
//            LinkageDiseqDisplayPlugin ldPlot = new LinkageDiseqDisplayPlugin(null, true);
//            ldPlot.performFunction(data);
//        }
    }

    public void actionPerformed(ActionEvent e) {
        enableButtons();
    }

    private void enableButtons() {
        if (panel3.isLDChecked() || panel3.isGLMChecked() || panel3.isMLMChecked() || panel3.isImputeChecked()) {
            getWizard().setNextFinishButtonEnabled(true);
        }
        else {
            getWizard().setNextFinishButtonEnabled(false);
        }
    }
}
