package net.maizegenetics.baseplugins;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;
import java.util.ArrayList;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;

import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.Trait;

public class FilterTraitsDialog extends JDialog implements ActionListener, TableModelListener {

    private Phenotype myPhenotype;
    private TraitTableModel traitModel;
    private JScrollPane jsp;
    private JTable traitTable;
    private boolean clickedOK = false;
    public final static String CMD_EXCLUDE = "exclude";
    public final static String CMD_INCLUDE = "include";
    public final static String CMD_EXCLUDE_ALL = "excludeall";
    public final static String CMD_INCLUDE_ALL = "includeall";
    public final static String CMD_OK = "ok";
    public final static String CMD_CANCEL = "cancel";

    public FilterTraitsDialog(Frame parent, Phenotype aPhenotype) {
        super(parent);
        myPhenotype = aPhenotype;
        traitModel = new TraitTableModel(myPhenotype);
        init();
    }

    private void init() {
        setTitle("Filter Traits / Modify Trait Properties");
        setSize(new Dimension(600, 600));
        setModal(true);
        Container contentPane = getContentPane();
        contentPane.setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();

        //add the table
        traitModel = new TraitTableModel(myPhenotype);
        traitTable = new JTable(traitModel);
        jsp = new JScrollPane(traitTable);

        //create a combo box editor for type
        JComboBox comboType = new JComboBox();
        comboType.addItem(Trait.TYPE_DATA);
        comboType.addItem(Trait.TYPE_COVARIATE);
        comboType.addItem(Trait.TYPE_FACTOR);
        comboType.addItem(Trait.TYPE_MARKER);
        int typeColNumber = traitModel.getTypeColumnNumber();

        TableColumn typeColumn = traitTable.getColumnModel().getColumn(typeColNumber);
        typeColumn.setCellEditor(new DefaultCellEditor(comboType));

        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridheight = 1;
        gbc.gridwidth = 2;
        contentPane.add(jsp, gbc);

        //add buttons
        JButton btnExclude = new JButton("Exclude Selected");
        btnExclude.setActionCommand(CMD_EXCLUDE);
        btnExclude.addActionListener(this);
        JButton btnInclude = new JButton("Include Selected");
        btnInclude.setActionCommand(CMD_INCLUDE);
        btnInclude.addActionListener(this);
        JButton btnExcludeAll = new JButton("Exclude All");
        btnExcludeAll.setActionCommand(CMD_EXCLUDE_ALL);
        btnExcludeAll.addActionListener(this);
        JButton btnIncludeAll = new JButton("Include All");
        btnIncludeAll.setActionCommand(CMD_INCLUDE_ALL);
        btnIncludeAll.addActionListener(this);
        JButton btnOK = new JButton("OK");
        btnOK.setActionCommand(CMD_OK);
        btnOK.addActionListener(this);
        JButton btnCancel = new JButton("Cancel");
        btnCancel.setActionCommand(CMD_CANCEL);
        btnCancel.addActionListener(this);

        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 1;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.weightx = 1;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        contentPane.add(btnExclude, gbc);

        gbc.gridx = 1;
        gbc.gridy = 1;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.WEST;
        contentPane.add(btnInclude, gbc);

        gbc.gridx = 0;
        gbc.gridy = 2;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.EAST;
        contentPane.add(btnExcludeAll, gbc);

        gbc.gridx = 1;
        gbc.gridy = 2;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.WEST;
        contentPane.add(btnIncludeAll, gbc);

        gbc.gridx = 0;
        gbc.gridy = 3;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.EAST;
        contentPane.add(btnOK, gbc);

        gbc.gridx = 1;
        gbc.gridy = 3;
        gbc.insets = new Insets(5, 5, 5, 5); //top, left, bottom, right
        gbc.anchor = GridBagConstraints.WEST;
        contentPane.add(btnCancel, gbc);
    }

    /**
     * @return an index of the traits to be included in the filter phenotype
     */
    public int[] getIncludedTraits() {
        if (!clickedOK) {
            return null;
        }
        Boolean[] isIncluded = traitModel.include;
        int howmany = 0;
        for (Boolean include : isIncluded) {
            if (include) {
                howmany++;
            }
        }
        int[] includedTrait = new int[howmany];
        int traitCount = 0;
        int includedCount = 0;
        for (Boolean include : isIncluded) {
            if (include) {
                includedTrait[includedCount++] = traitCount;
            }
            traitCount++;
        }
        return includedTrait;
    }

    /**
     * @return trait types for all the traits in the original phenotype, not just the included ones
     */
    public String[] getTraitTypes() {
        int n = traitModel.traitList.size();
        String[] types = new String[n];
        for (int i = 0; i < n; i++) {
            types[i] = traitModel.traitList.get(i).getType();
        }
        return types;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equals(CMD_EXCLUDE)) {
            traitModel.excludeSome(traitTable.getSelectedRows());
        } else if (e.getActionCommand().equals(CMD_INCLUDE)) {
            traitModel.includeSome(traitTable.getSelectedRows());
        } else if (e.getActionCommand().equals(CMD_EXCLUDE_ALL)) {
            traitModel.excludeAll();
        } else if (e.getActionCommand().equals(CMD_INCLUDE_ALL)) {
            traitModel.includeAll();
        } else if (e.getActionCommand().equals(CMD_OK)) {
            setVisible(false);
            clickedOK = true;
        } else if (e.getActionCommand().equals(CMD_CANCEL)) {
            setVisible(false);
            clickedOK = false;
        }
    }

    @Override
    public void tableChanged(TableModelEvent e) {
        System.out.println("table model changed");
    }

    public boolean getClickedOK() {
        return clickedOK;
    }
}

class TraitTableModel extends AbstractTableModel {

    ArrayList<Trait> traitList;
    Phenotype myPhenotype;
    Boolean[] include;
    String[] colName;
    int numberOfTraits;
    int numberOfFactors;
    int numberOfColumns;
    int typeColumnNumber;

    TraitTableModel(Phenotype aPhenotype) {
        super();
        myPhenotype = aPhenotype;
        traitList = new ArrayList<Trait>();
        List<Trait> oldList = myPhenotype.getTraits();
        for (Trait trait : oldList) {
            traitList.add(Trait.getInstance(trait));
        }
        numberOfTraits = traitList.size();
        numberOfFactors = myPhenotype.getNumberOfFactors();
        numberOfColumns = 4 + numberOfFactors;
        setColumnNames();
        include = new Boolean[numberOfTraits];
        for (int i = 0; i < numberOfTraits; i++) {
            include[i] = true;
        }
    }

    private void setColumnNames() {
        colName = new String[numberOfColumns];
        int col = 0;
        colName[col++] = "Trait";
        for (String factorname : myPhenotype.getFactorNameCopy()) {
            colName[col++] = factorname;
        }
        typeColumnNumber = col;
        colName[col++] = "Type";
        colName[col++] = "Discrete";
        colName[col++] = "Include";
    }

    @Override
    public int getColumnCount() {
        return colName.length;
    }

    @Override
    public int getRowCount() {
        return include.length;
    }

    @Override
    public Object getValueAt(int row, int col) {
        if (col == 3 + numberOfFactors) {
            return include[row];
        }
        if (col == 0) {
            return traitList.get(row).getName();
        }
        if (col == 1 + numberOfFactors) {
            return traitList.get(row).getType();
        }
        if (col == 2 + numberOfFactors) {
            return traitList.get(row).isDiscrete();
        } else {
            String factorName = myPhenotype.getFactorName(col - 1);
            Object factorValue = myPhenotype.getTrait(row).getProperty(factorName);
            if (factorValue == null) {
                return "";
            }
            return factorValue;
        }
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        if (columnIndex < 2 + numberOfFactors) {
            return String.class;
        } else {
            return Boolean.class;
        }
    }

    @Override
    public String getColumnName(int column) {
        return colName[column];
    }

    @Override
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        if (columnIndex == typeColumnNumber || columnIndex == typeColumnNumber + 2) {
            return true;
        }
        return false;
    }

    @Override
    public void setValueAt(Object value, int rowIndex, int columnIndex) {
        if (columnIndex == 3 + numberOfFactors) {
            include[rowIndex] = (Boolean) value;
        }
        if (columnIndex == 1 + numberOfFactors) {
            traitList.get(rowIndex).setType(value.toString());
        }
    }

    public void excludeAll() {
        for (int i = 0; i < numberOfTraits; i++) {
            include[i] = Boolean.FALSE;
        }
        fireTableDataChanged();
    }

    public void includeAll() {
        for (int i = 0; i < numberOfTraits; i++) {
            include[i] = Boolean.TRUE;
        }
        fireTableDataChanged();
    }

    public void excludeSome(int[] index) {
        for (int i : index) {
            include[i] = Boolean.FALSE;
        }
        fireTableDataChanged();
    }

    public void includeSome(int[] index) {
        for (int i : index) {
            include[i] = Boolean.TRUE;
        }
        fireTableDataChanged();
    }

    public int getTypeColumnNumber() {
        return typeColumnNumber;
    }
}
