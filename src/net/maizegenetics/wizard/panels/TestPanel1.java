package net.maizegenetics.wizard.panels;

import java.awt.*;
import java.net.*;
import javax.swing.*;
import javax.swing.border.*;


public class TestPanel1 extends JPanel {
 
    private JLabel blankSpace;
    private JLabel jLabel1;
    private JLabel jLabel2;
    private JLabel jLabel3;
    private JLabel jLabel4;
    private JLabel jLabel5;
    private JLabel jLabel6;
    private JLabel jLabel7;
    private JLabel jLabel8;
    private JLabel jLabel9;

    private JLabel welcomeTitle;
    private JPanel contentPanel;
    
    //private JLabel iconLabel;
    //private ImageIcon icon;
    
    public TestPanel1() {
        
        //iconLabel = new JLabel();
        contentPanel = getContentPanel();
        //contentPanel.setMinimumSize(new Dimension(640, 480));
        contentPanel.setBorder(new EmptyBorder(new Insets(10, 10, 10, 10)));

        //icon = getImageIcon();

        setLayout(new java.awt.BorderLayout());

        //if (icon != null)
            //iconLabel.setIcon(icon);
        
        //iconLabel.setBorder(new EtchedBorder(EtchedBorder.RAISED));
        
        //add(iconLabel, BorderLayout.WEST);
        
        JPanel secondaryPanel = new JPanel();
        secondaryPanel.add(contentPanel, BorderLayout.NORTH);
        add(secondaryPanel, BorderLayout.CENTER);
    }
    
    
    private JPanel getContentPanel() {
        
        JPanel contentPanel1 = new JPanel();
        JPanel jPanel1 = new JPanel();
        
        welcomeTitle = new JLabel();
        blankSpace = new JLabel();
        jLabel1 = new JLabel();
        jLabel2 = new JLabel();
        jLabel3 = new JLabel();
        jLabel4 = new JLabel();
        jLabel5 = new JLabel();
        jLabel6 = new JLabel();
        jLabel7 = new JLabel();
        jLabel8 = new JLabel();
        jLabel9 = new JLabel();

        contentPanel1.setLayout(new java.awt.BorderLayout());

        welcomeTitle.setFont(new java.awt.Font("MS Sans Serif", Font.BOLD, 14));
        welcomeTitle.setText("Welcome to the TASSEL Wizard!");
        contentPanel1.add(welcomeTitle, java.awt.BorderLayout.NORTH);

        jPanel1.setLayout(new java.awt.GridLayout(0, 1));

        jPanel1.add(blankSpace);
        //jLabel1.setText("This is an example of a wizard dialog, which allows a user to traverse");
        jLabel1.setText("This Wizard is designed to help new TASSEL users to quickly and easily");
        jPanel1.add(jLabel1);
        //jLabel2.setText("a number of panels (while entering data) until the wizard has enough ");
        jLabel2.setText("begin using tools available for data analysis without having to first");
        jPanel1.add(jLabel2);
        //jLabel3.setText("information to perform whatever end function is necessary. Note that");
        jLabel3.setText("become familiar with the main TASSEL user interface. If you would like");
        jPanel1.add(jLabel3);
        //jLabel4.setText("panels are not necessarily ordered in a linear fashion, but instead in");
        jLabel4.setText("to proceed directly to the TASSEL user interface, hit \"Cancel\" now.");
        jPanel1.add(jLabel4);
        //jLabel5.setText("a tree-like manner (e.g., there may be more than one panel with a");
        jPanel1.add(jLabel5);
        //jLabel7.setText("Otherwise plase hit \"Next\" to continue the TASSEL Wizard....");
        jPanel1.add(jLabel6);
        //jLabel6.setText("traverse the path). That's not the case with this example, however.");
        jPanel1.add(jLabel7);
        jLabel8.setText("Otherwise plase hit \"Next\" to continue the TASSEL Wizard....");
        jPanel1.add(jLabel8);
        //jLabel9.setText("Press the 'Next' button to continue....");
        jPanel1.add(jLabel9);

        contentPanel1.add(jPanel1, java.awt.BorderLayout.CENTER);

        return contentPanel1;
        
    }

//    private ImageIcon getImageIcon() {
//        return new ImageIcon((URL)getResource("clouds.jpg"));
//    }
    
    private Object getResource(String key) {

        URL url = null;
        String name = key;

        if (name != null) {

            try {
                Class c = Class.forName("net.maizegenetics.wizard.panels.Main");
                url = c.getResource(name);
            } catch (ClassNotFoundException cnfe) {
                System.err.println("Unable to find Main class");
            }
            return url;
        } else
            return null;

    }
 
}
