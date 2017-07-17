package net.maizegenetics.baseplugins.numericaltransform;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;

/**
 */
public class NumericalGenotypePanel extends JPanel {

    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JRadioButton myCollapse = new JRadioButton();
    private JRadioButton mySeparate = new JRadioButton();

    public NumericalGenotypePanel(Datum theDatum) throws Exception {

        try {
            if (!(theDatum.getData() instanceof Alignment)) {
                throw new Exception("Must be a Genotype");
            }

            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private void jbInit() throws Exception {

        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));

        //Radio Buttons
        myCollapse.setText("Collapse Non Major Alleles");
        mySeparate.setText("Separate Alleles");

        myButtonGroup.add(myCollapse);
        myButtonGroup.add(mySeparate);
        myButtonGroup.setSelected(myCollapse.getModel(), true);

        panel.add(myCollapse);
        panel.add(mySeparate);

        this.add(panel);

    }

    public boolean isCollapseSelected() {
        return myCollapse.isSelected();
    }

    public boolean isSeparateSelected() {
        return mySeparate.isSelected();
    }
}


