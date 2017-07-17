/*
 * PluginProgressUnit
 */
package net.maizegenetics.progress;

import java.awt.FlowLayout;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.net.URL;

import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginListener;
import net.maizegenetics.plugindef.ThreadedPluginListener;

/**
 *
 * @author terry
 */
public class PluginProgressUnit extends JPanel implements PluginListener {

    private final int LABEL_WIDTH = 20;
    private final JProgressBar myProgress;
    private final JButton myCancelButton;
    private final Plugin myPlugin;
    private final Plugin myCancelPlugin;
    private int myCurrentValue = 0;

    public PluginProgressUnit(Plugin plugin, int level) {
        this(plugin, level, true);
    }

    public PluginProgressUnit(Plugin plugin, int level, boolean cancelButton) {
        this(plugin, level, cancelButton, null);
    }

    public PluginProgressUnit(Plugin plugin, int level, boolean cancelButton, Plugin cancelPlugin) {

        myPlugin = plugin;

        if (cancelPlugin == null) {
            myCancelPlugin = myPlugin;
        } else {
            myCancelPlugin = cancelPlugin;
        }

        BoxLayout layout = new BoxLayout(this, BoxLayout.Y_AXIS);
        setLayout(layout);
        setAlignmentX(JPanel.TOP_ALIGNMENT);

        JPanel top = new JPanel();
        top.setLayout(new FlowLayout(FlowLayout.LEFT));

        myProgress = new JProgressBar(0, 100);
        setProgress(0);

        if (cancelButton) {

            URL imageURL = PluginProgressUnit.class.getResource("cancel.gif");
            if (imageURL == null) {
                myCancelButton = new JButton("Cancel");
            } else {
                myCancelButton = new JButton(new ImageIcon(imageURL));
            }
            myCancelButton.setToolTipText("Cancel");

            myCancelButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    if (!myPlugin.cancel()) {
                        stopThread(myCancelPlugin);
                    }
                    myCancelButton.setEnabled(false);
                    setProgress(100);
                }
            });
        } else {
            myCancelButton = null;
        }

        StringBuilder builder = new StringBuilder();
        for (int i = 0, n = level * 4; i < n; i++) {
            builder.append(" ");
        }
        String buttonName = myPlugin.getButtonName();
        if (buttonName.length() <= LABEL_WIDTH) {
            builder.append(buttonName);
            for (int i = 0, n = LABEL_WIDTH - buttonName.length(); i < n; i++) {
                builder.append(" ");
            }
        } else {
            builder.append(buttonName.substring(0, LABEL_WIDTH));
        }

        JLabel pluginLabel = new JLabel(builder.toString());
        pluginLabel.setFont(new Font("Monospaced", Font.PLAIN, 14));

        top.add(pluginLabel);
        top.add(myProgress);
        if (myCancelButton != null) {
            top.add(myCancelButton);
        }

        add(top);

        myPlugin.addListener(this);

    }

    private static void stopThread(Plugin plugin) {
        ThreadGroup root = Thread.currentThread().getThreadGroup().getParent();
        while (root.getParent() != null) {
            root = root.getParent();
        }

        visit(root, 0, plugin);
    }

    private static void visit(ThreadGroup group, int level, Plugin plugin) {
        // Get threads in `group'
        int numThreads = group.activeCount();
        Thread[] threads = new Thread[numThreads * 2];
        numThreads = group.enumerate(threads, false);

        // Enumerate each thread in `group'
        for (int i = 0; i < numThreads; i++) {
            // Get thread
            Thread thread = threads[i];
            if (thread instanceof ThreadedPluginListener) {
                if (plugin == ((ThreadedPluginListener) thread).getPluginListener()) {
                    try {
                        thread.stop();
                    } catch (Exception e) {
                        // do nothing
                    }
                    ProgressPanel.getInstance().cleanProgressUnit(plugin);
                }
            }
        }

        // Get thread subgroups of `group'
        int numGroups = group.activeGroupCount();
        ThreadGroup[] groups = new ThreadGroup[numGroups * 2];
        numGroups = group.enumerate(groups, false);

        // Recursively visit each subgroup
        for (int i = 0; i < numGroups; i++) {
            visit(groups[i], level + 1, plugin);
        }
    }

    public void dataSetReturned(PluginEvent event) {
        setProgress(100);
    }

    public void progress(PluginEvent event) {
        DataSet dataSet = (DataSet) event.getSource();
        List<Datum> datum = dataSet.getDataOfType(Integer.class);
        if (datum.size() != 1) {
            throw new IllegalArgumentException("PluginProgressUnit: progress: should be only one interger percent value.");
        }
        Integer percent = (Integer) ((Datum) datum.get(0)).getData();
        setProgress(percent.intValue());
    }

    private void setProgress(int num) {
        if (myCurrentValue == num) {
            return;
        }
        myCurrentValue = num;
        myProgress.setValue(num);
        if (num == 100) {
            ProgressPanel.getInstance().removeProgressUnit(this);
        }
    }
}
