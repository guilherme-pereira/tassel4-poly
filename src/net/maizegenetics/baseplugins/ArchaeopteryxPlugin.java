package net.maizegenetics.baseplugins;

import java.awt.Frame;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;
import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

import net.maizegenetics.pal.tree.Node;
import net.maizegenetics.pal.tree.Tree;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class ArchaeopteryxPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ArchaeopteryxPlugin.class);

    public ArchaeopteryxPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        try {
            List<Datum> treeInList = input.getDataOfType(Tree.class);
            if (treeInList.size() != 1) {
                String message = "Invalid selection.  Please select one tree.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }

            Tree myTree = (Tree) treeInList.get(0).getData();

            makeTree(myTree, isInteractive(), getParentFrame());
            return null;
        } finally {
            fireProgress(100);
        }
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = TreeDisplayPlugin.class.getResource("images/Tree.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Archaeopteryx Tree";
    }

    @Override
    public String getToolTipText() {
        return "Archaeopteryx Tree Viewer";
    }

    public static void makeTree(Tree theTree, boolean isInteractive, Frame frame) {

        class nodePair {

            PhylogenyNode fNode;
            Node pNode;

            nodePair(PhylogenyNode node1, Node node2) {
                fNode = node1;
                pNode = node2;
            }
        }

        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        Node palroot = theTree.getRoot();
        root.setName(palroot.getIdentifier().getName());
        LinkedList<nodePair> nodeList = new LinkedList<nodePair>();
        nodeList.add(new nodePair(root, palroot));

        while (!nodeList.isEmpty()) {
            nodePair pair = nodeList.remove();
            int n = pair.pNode.getChildCount();
            for (int i = 0; i < n; i++) {
                PhylogenyNode newNode = new PhylogenyNode();
                pair.fNode.addAsChild(newNode);
                newNode.setName(pair.pNode.getChild(i).getIdentifier().getName());
                newNode.setDistanceToParent(pair.pNode.getChild(i).getBranchLength());
                nodeList.add(new nodePair(newNode, pair.pNode.getChild(i)));
            }
        }
        phy.setRoot(root);
        phy.setRooted(true);

        if (!isInteractive) {
            myLogger.info("Archaeopteryx plugin non-interactive mode is not implemented.");
            return;
        }

        Archaeopteryx.createApplication(phy).setLocationRelativeTo(frame);
    }
}
