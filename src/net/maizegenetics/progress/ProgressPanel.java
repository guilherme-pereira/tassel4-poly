/*
 * ProgressPanel
 */
package net.maizegenetics.progress;

import java.awt.BorderLayout;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.Plugin;

/**
 *
 * @author terry
 */
public class ProgressPanel extends JPanel {

    private static ProgressPanel SINGLETON = null;
    private final JPanel myMainPane;
    private final Map myPlugins = new HashMap();
    private final Map myCancelPlugins = new HashMap();

    private ProgressPanel() {

        setLayout(new BorderLayout());

        myMainPane = new JPanel();
        BoxLayout layout = new BoxLayout(myMainPane, BoxLayout.Y_AXIS);
        myMainPane.setLayout(layout);
        myMainPane.setAlignmentX(JPanel.TOP_ALIGNMENT);

        JScrollPane scrollPane = new JScrollPane();
        scrollPane.getViewport().add(myMainPane);

        add(scrollPane, BorderLayout.CENTER);

    }

    public static ProgressPanel getInstance() {
        if (SINGLETON == null) {
            SINGLETON = new ProgressPanel();
        }
        return SINGLETON;
    }

    public void addPipelineSegment(List<Plugin> plugins) {

        if ((plugins == null) || (plugins.size() == 0)) {
            return;
        }

        for (int i = 0; i < plugins.size(); i++) {
            addPlugin(plugins.get(i), true, plugins.get(0));
        }

    }

    public void addPlugin(Plugin plugin) {
        addPlugin(plugin, true, null);
    }

    public void addPlugin(Plugin plugin, boolean cancelButton, Plugin cancelPlugin) {

        Integer level = new Integer(0);
        Iterator itr = ((AbstractPlugin) plugin).getInputs().iterator();
        while (itr.hasNext()) {
            Plugin current = (Plugin) itr.next();
            Integer inputLevel = (Integer) myPlugins.get(current);
            if ((inputLevel != null) && (inputLevel.compareTo(level) >= 0)) {
                level = inputLevel + 1;
            }
        }

        myPlugins.put(plugin, level);

        PluginProgressUnit current = new PluginProgressUnit(plugin, level.intValue(), cancelButton, cancelPlugin);

        if (cancelPlugin != null) {
            myCancelPlugins.put(current, cancelPlugin);
        }

        myMainPane.add(current);

        revalidate();

    }

    public void removeProgressUnit(PluginProgressUnit unit) {
        synchronized (myCancelPlugins) {
            Plugin cancelPlugin = (Plugin) myCancelPlugins.get(unit);
            if (cancelPlugin != null) {
                myCancelPlugins.remove(unit);
            }
            myMainPane.remove(unit);
        }
        revalidate();
        repaint();
    }

    public void cleanProgressUnit(Plugin plugin) {
        synchronized (myCancelPlugins) {
            Iterator itr = myCancelPlugins.entrySet().iterator();
            while (itr.hasNext()) {
                Map.Entry current = (Map.Entry) itr.next();
                Plugin tempPlugin = (Plugin) current.getValue();
                if (tempPlugin == plugin) {
                    itr.remove();
                    myMainPane.remove((PluginProgressUnit) current.getKey());
                }
            }
        }
        revalidate();
        repaint();
    }
}
