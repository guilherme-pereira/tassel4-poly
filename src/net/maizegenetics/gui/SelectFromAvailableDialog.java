/*
 * SelectFromAvailableDialog
 */
package net.maizegenetics.gui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;

import java.util.Arrays;

import javax.swing.AbstractAction;
import javax.swing.AbstractListModel;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JRootPane;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.ScrollPaneConstants;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author terry
 */
public class SelectFromAvailableDialog extends JDialog {

    private final static int LIST_ROW_COUNT = 20;
    private JList myAvailable;
    private JList mySelected;
    private AbstractAvailableListModel myAvailableListModel;
    private SelectedListModel mySelectedListModel;
    private boolean myIsCanceled = false;
    private boolean myIsCaptureSelected = true;

    public SelectFromAvailableDialog(Frame frame, String title, AbstractAvailableListModel availableModel) {
        super(frame, true);

        myAvailableListModel = availableModel;
        myAvailable = new JList(myAvailableListModel);
        myAvailable.setPrototypeCellValue("XXXXXXXXXXXXXXXXXXXX");

        mySelectedListModel = new SelectedListModel();
        mySelected = new JList(mySelectedListModel);
        mySelected.setPrototypeCellValue("XXXXXXXXXXXXXXXXXXXX");

        setTitle(title);
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);

        Container contentPane = getContentPane();

        BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
        contentPane.setLayout(layout);

        JPanel whole = new JPanel(new BorderLayout());

        whole.add(getMain(), BorderLayout.CENTER);
        whole.add(getBottomPanel(), BorderLayout.SOUTH);

        contentPane.add(whole);

        pack();
        validate();

        setResizable(true);

    }

    private JPanel getMain() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.X_AXIS);
        result.setLayout(layout);

        result.add(Box.createRigidArea(new Dimension(15, 1)));

        JPanel available = getLabeledListWSearch("Available", myAvailable);
        result.add(available);

        result.add(Box.createRigidArea(new Dimension(15, 1)));
        result.add(getAddButton());
        result.add(Box.createRigidArea(new Dimension(15, 1)));

        JPanel selected = getLabeledList("Selected", mySelected);
        result.add(selected);

        result.add(Box.createRigidArea(new Dimension(15, 1)));

        return result;

    }

    private JPanel getLabeledList(String label, JList list) {
        JPanel result = new JPanel();
        result.setLayout(new BorderLayout());
        result.setPreferredSize(new Dimension(150, 200));
        result.add(new JLabel(label), BorderLayout.NORTH);
        list.setVisibleRowCount(LIST_ROW_COUNT);
        list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        JScrollPane scrollPane = new JScrollPane(list);
        scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        result.add(scrollPane, BorderLayout.CENTER);
        return result;
    }

    private JPanel getLabeledListWSearch(String label, JList list) {
        JPanel result = new JPanel();
        result.setLayout(new BorderLayout());
        result.setPreferredSize(new Dimension(150, 200));
        result.add(new JLabel(label), BorderLayout.NORTH);
        JTextField search = new JTextField();
        CaretListener caret = new CaretListener() {
            @Override
            public void caretUpdate(CaretEvent e) {
                String text = ((JTextField) e.getSource()).getText();
                myAvailableListModel.setShown(text);
            }
        };
        search.addCaretListener(caret);
        result.add(search, BorderLayout.SOUTH);
        list.setVisibleRowCount(LIST_ROW_COUNT);
        list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        JScrollPane scrollPane = new JScrollPane(list);
        scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        result.add(scrollPane, BorderLayout.CENTER);
        return result;
    }

    private JButton getAddButton() {
        JButton result = new JButton("Add ->");
        result.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                int[] selected = myAvailable.getSelectedIndices();
                selected = myAvailableListModel.translateToRealIndices(selected);
                for (int i = 0, n = selected.length; i < n; i++) {
                    mySelectedListModel.add(selected[i]);
                }
            }
        });
        return result;
    }

    private JPanel getBottomPanel() {
        JPanel result = new JPanel(new FlowLayout());

        result.add(getSelectedButton());
        result.add(Box.createRigidArea(new Dimension(15, 1)));
        result.add(getUnselectedButton());
        result.add(Box.createRigidArea(new Dimension(15, 1)));
        result.add(getRemoveButton());
        result.add(Box.createRigidArea(new Dimension(15, 1)));
        result.add(getCancelButton());

        return result;
    }

    private JButton getSelectedButton() {
        JButton okButton = new JButton("Capture Selected");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                setVisible(false);
                myIsCanceled = false;
                myIsCaptureSelected = true;
            }
        });
        return okButton;
    }

    private JButton getUnselectedButton() {
        JButton okButton = new JButton("Capture Unselected");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                setVisible(false);
                myIsCanceled = false;
                myIsCaptureSelected = false;
            }
        });
        return okButton;
    }

    private JButton getCancelButton() {
        JButton okButton = new JButton("Cancel");
        okButton.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                setVisible(false);
                myIsCanceled = true;
            }
        });
        return okButton;
    }

    private JButton getRemoveButton() {
        JButton result = new JButton("Remove");
        result.addActionListener(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                int[] indices = mySelected.getSelectedIndices();
                Arrays.sort(indices);
                for (int i = indices.length - 1; i >= 0; i--) {
                    mySelectedListModel.remove(indices[i]);
                    mySelected.clearSelection();
                }
            }
        });
        return result;
    }

    public boolean isCanceled() {
        return myIsCanceled;
    }

    public int[] getDesiredIndices() {

        if (myIsCaptureSelected) {
            return mySelectedListModel.getBitSet().getIndicesOfSetBits();
        } else {
            BitSet temp = mySelectedListModel.getBitSet();
            temp.flip(0, myAvailableListModel.getRealSize());
            return temp.getIndicesOfSetBits();
        }

    }

    public class SelectedListModel extends AbstractListModel {

        private BitSet mySelectedIndices;

        public SelectedListModel() {
            mySelectedIndices = new OpenBitSet(myAvailableListModel.getRealSize());
        }

        @Override
        public int getSize() {
            return (int) mySelectedIndices.cardinality();
        }

        @Override
        public Object getElementAt(int index) {
            return myAvailableListModel.getRealElementAt(mySelectedIndices.indexOfNthSetBit(index + 1));
        }

        public void add(int index) {
            mySelectedIndices.fastSet(index);
            fireContentsChanged(this, 0, getSize() - 1);
        }

        public void remove(int index) {
            mySelectedIndices.fastClear(mySelectedIndices.indexOfNthSetBit(index + 1));
            fireContentsChanged(this, 0, getSize() - 1);
        }

        public BitSet getBitSet() {
            return mySelectedIndices;
        }
    }
}
