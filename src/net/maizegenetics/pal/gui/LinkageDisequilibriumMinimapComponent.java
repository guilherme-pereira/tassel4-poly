/**
 * Created 2/24/12
 *
 */

package net.maizegenetics.pal.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import javax.swing.JComponent;

import net.maizegenetics.pal.popgen.LinkageDisequilibrium;

/**
 *
 * @author yz79
 */
public class LinkageDisequilibriumMinimapComponent extends JComponent {

    //LD attributes
    LinkageDisequilibrium myLD;
    int myNumSites;

    //Viewable window attributes
    int myWindowSize;
    int myWindowX;
    int myWindowY;

    //Minimap attributes
    final int myMapSize = 100;
    final int myDiagonalWidth = 3; //pixels on each side counting center diagonal

    //Current viewable window with respect to minimap, on a scale of 0 - 1
    double myXStart;
    double myYStart;

    double myXEnd;
    double myYEnd;

    //Current viewable window in pixels
    int myXStartPos;
    int myYStartPos;

    int myXEndPos;
    int myYEndPos;


    public LinkageDisequilibriumMinimapComponent(LinkageDisequilibrium theLD, int windowSize, int windowX, int windowY) {

        myLD = theLD;
        myNumSites = theLD.getSiteCount();
        myWindowSize = windowSize;
        myWindowX = windowX;
        myWindowY = windowY;

        calculateStartEndPositions();
        calculateCoordinates();

        try {
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * resets myWindowSize
     * @param windowSize
     * @param windowX
     * @param windowY
     */
    public void setWindowSize(int windowSize, int windowX, int windowY) {

        myWindowSize = windowSize;
        myWindowX = windowX;
        myWindowY = windowY;

        calculateStartEndPositions();
        calculateCoordinates();
    }

    /**
     * resets myWindowX
     * @param windowX
     */
    public void setWindowX(int windowX) {

        myWindowX = windowX;

        calculateStartEndPositions();
        calculateCoordinates();
    }

    /**
     * resets myWindowY
     * @param windowY
     */
    public void setWindowY(int windowY) {

        myWindowY = windowY;

        calculateStartEndPositions();
        calculateCoordinates();
    }

    /**
     * Caluclates start and end positions of viewable window on a scale of 0 - 1
     */
    private void calculateStartEndPositions() {

        myXStart = Math.floor(myWindowX-myWindowSize/2.0)/myNumSites;
        myYStart = Math.floor(myWindowY-myWindowSize/2.0)/myNumSites;

        myXEnd = (Math.floor(myWindowX-myWindowSize/2.0)+myWindowSize)/myNumSites;
        myYEnd = (Math.floor(myWindowY-myWindowSize/2.0)+myWindowSize)/myNumSites;
    }

    /**
     * Calculates start and end positions of viewable window in physical position
     * compensates for box becoming too small, minimum size = 2x2
     */
    private void calculateCoordinates() {

        myXStartPos = (int)(myXStart*myMapSize);
        myYStartPos = (int)(myYStart*myMapSize);

        myXEndPos = (int)(myXEnd*myMapSize)-1;
        myYEndPos = (int)(myYEnd*myMapSize)-1;

        if (!(myXStartPos < myMapSize-1)) {
            myXStartPos = myMapSize - 2;
            myXEndPos = myMapSize - 1;
        }

        if (!(myYStartPos < myMapSize-1)) {
            myYStartPos = myMapSize - 2;
            myYEndPos = myMapSize - 1;
        }

        if (myXEndPos <= myXStartPos) {
            myXEndPos = myXStartPos + 1;
        }

        if (myYEndPos <= myYStartPos) {
            myYEndPos = myYStartPos + 1;
        }
    }

    private void jbInit() throws Exception {

        setBackground(Color.white);
        setPreferredSize(new Dimension(myMapSize, myMapSize));
    }

    protected void paintComponent(Graphics g) {

        Dimension d = getSize();

        g.setColor(Color.white);
        g.fillRect(0, 0, myMapSize, myMapSize);
        paintDiagonal(g);

        //paints red box
        g.setColor(Color.red);
        g.drawRect(myXStartPos, myYStartPos, myXEndPos-myXStartPos, myYEndPos-myYStartPos);
    }

    public void paint(Graphics g) {
        paintComponent(g);
    }

    /**
     * Adds diagonal
     * @param g
     */
    private void paintDiagonal(Graphics g) {
        g.setColor(Color.black);
        for (int i = 0; i < myDiagonalWidth; i++) {
            g.drawLine(0, 0+i, myMapSize-1-i, myMapSize-1);
            g.drawLine(0+i, 0, myMapSize-1, myMapSize-1-i);
        }
    }

    /**
     * Returns map size
     * @return
     */
    public int getMapSize() {
        return myMapSize;
    }
}
