package net.maizegenetics.wizard.panels;

import net.maizegenetics.wizard.WizardPanelDescriptor;
//import net.maizegenetics.wizard.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import net.maizegenetics.baseplugins.LinkageDiseqDisplayPlugin;


public class TestPanel4Descriptor extends WizardPanelDescriptor implements ActionListener {

    public static final String IDENTIFIER = "SERVER_CONNECT_PANEL";

    TestPanel4 panel4;

    public TestPanel4Descriptor() {

        panel4 = new TestPanel4();
        setPanelDescriptorIdentifier(IDENTIFIER);
        setPanelComponent(panel4);

    }

    public Object getNextPanelDescriptor() {
        return FINISH;
    }

    public Object getBackPanelDescriptor() {
        return TestPanel3Descriptor.IDENTIFIER;
    }


    public void aboutToDisplayPanel() {

//        panel4.setProgressValue(0);
//        panel4.setProgressText("Connecting to Server...");

        //getWizard().setNextFinishButtonEnabled(false);
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
//                    panel4.setProgressValue(25);
//                    panel4.setProgressText("Server Connection Established");
//                    Thread.sleep(500);
//                    panel4.setProgressValue(50);
//                    panel4.setProgressText("Transmitting Data...");
//                    Thread.sleep(3000);
//                    panel4.setProgressValue(75);
//                    panel4.setProgressText("Receiving Acknowledgement...");
//                    Thread.sleep(1000);
//                    panel4.setProgressValue(100);
//                    panel4.setProgressText("Data Successfully Transmitted");
//
//                    getWizard().setNextFinishButtonEnabled(true);
//                    getWizard().setBackButtonEnabled(true);
//
//                } catch (InterruptedException e) {
//
//                    panel4.setProgressValue(0);
//                    panel4.setProgressText("An Error Has Occurred");
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
        if (panel4.isLDPlotChecked()) {
            LinkageDiseqDisplayPlugin ldPlot = new LinkageDiseqDisplayPlugin(null, true);
            ldPlot.performFunction(data);
        }
    }

    public void actionPerformed(ActionEvent e) {
//        enableButtons();
    }

//    private void enableButtons() {
//        if (panel3.isLDChecked() || panel3.isGLMChecked() || panel3.isMLMChecked() || panel3.isImputeChecked()) {
//            getWizard().setNextFinishButtonEnabled(true);
//        }
//        else {
//            getWizard().setNextFinishButtonEnabled(false);
//        }
//    }

}
