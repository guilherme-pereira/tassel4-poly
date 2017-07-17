/*
 * TasselPipelineXMLUtil
 */
package net.maizegenetics.pipeline;

import java.awt.Frame;

import java.io.File;
import java.io.IOException;

import java.lang.reflect.Constructor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import net.maizegenetics.util.Utils;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 *
 * @author terry
 */
public class TasselPipelineXMLUtil {

    /**
     * These are flags that create plugins.
     */
    private static final String[] TAG_STRINGS = new String[]{"-t", "-s",
        "-p", "-a", "-k", "-q", "-h", "-b", "-g", "-r", "-plink", "-flapjack", "-fasta", "-geneticMap",
        "-taxaJoinStrict", "-union", "-intersect", "-excludeLastTrait", "-mlm", "-glm", "-wd_csv", "-wd_tab",
        "-td_csv", "-td_tab", "-td_gui", "-wxls", "-ld", "-ldd", "-ck", "-gs", "-export",
        "-impute", "-filterAlign", "-includeTaxa", "-includeTaxaInFile", "-excludeTaxa"};
    private static final List TAG_FLAGS = new ArrayList(Arrays.asList(TAG_STRINGS));

    private TasselPipelineXMLUtil() {
        // Utility Class
    }

    public static void writeArgsAsXML(String filename, String[] args) {

        try {
            DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder docBuilder = docFactory.newDocumentBuilder();

            Document doc = docBuilder.newDocument();
            Element rootElement = doc.createElement("TasselPipeline");
            doc.appendChild(rootElement);
            createXML(doc, rootElement, args);

            TransformerFactory transformerFactory = TransformerFactory.newInstance();
            Transformer transformer = transformerFactory.newTransformer();
            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");
            transformer.setOutputProperty(OutputKeys.METHOD, "xml");
            DOMSource source = new DOMSource(doc);
            StreamResult result = new StreamResult(new File(filename));

            transformer.transform(source, result);
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
        }

    }

    private static void createXML(Document doc, Element element, String[] args) throws IOException {

        int index = 0;

        while (index < args.length) {
            String current = args[index];
            if (!isFork(current)) {
                throw new IllegalArgumentException("TasselPipelineXMLUtil: createXML: this flag should be either -fork, -combine, or -runfork: " + current);
            }
            Element newElement = createTag(doc, element, current);

            /*
            if (current.startsWith("-fork")) {
                String name = current.replaceFirst("-fork", "");
                newElement = createTag(doc, element, "-fork");
                newElement.setAttribute("name", name);
            } else if (current.startsWith("-runfork")) {
                String name = current.replaceFirst("-runfork", "");
                newElement = createTag(doc, element, "-runfork");
                newElement.setAttribute("name", name);
            } else if (current.startsWith("-combine")) {
                String name = current.replaceFirst("-combine", "");
                newElement = createTag(doc, element, "-combine");
                newElement.setAttribute("name", name);
            }
             */

            while (true) {
                index++;
                if (index >= args.length) {
                    break;
                }
                if (isFork(args[index])) {
                    break;
                } else if (isSelfDescribingPlugin(args[index])) {
                    index = createSelfDescribingPluginXML(doc, newElement, args, index);
                } else if (isModifier(args[index])) {
                    index = createString(doc, newElement, args, index);
                } else {
                    index = createXML(doc, newElement, args, index);
                }
            }

        }

    }

    private static int createXML(Document doc, Element element, String[] args, int index) throws IOException {

        String current = args[index];
        Element newElement = createTag(doc, element, current);

        while (true) {
            index++;
            if (index >= args.length) {
                break;
            }
            if (isModifier(args[index])) {
                index = createString(doc, newElement, args, index);
            } else {
                break;
            }
        }

        return index - 1;

    }

    private static int createSelfDescribingPluginXML(Document doc, Element element, String[] args, int index) throws IOException {

        String current = args[index];
        Element newElement = createTag(doc, element, current);

        while (true) {
            index++;
            if (index >= args.length) {
                break;
            }
            if (args[index].equalsIgnoreCase("-endPlugin")) {
                index++;
                break;
            }
            index = createString(doc, newElement, args, index);
        }

        return index - 1;

    }

    private static boolean isFork(String str) {

        if ((str.startsWith("-fork")) || (str.startsWith("-runfork")) || (str.startsWith("-combine"))) {
            return true;
        } else {
            return false;
        }

    }

    private static boolean isModifier(String str) {

        if (str.startsWith("-")) {

            if ((str.startsWith("-fork")) || (str.startsWith("-runfork")) || (str.startsWith("-combine"))) {
                return false;
            } else if (TAG_FLAGS.contains(str)) {
                return false;
            } else {
                return !isSelfDescribingPlugin(str);
            }

        } else {
            return true;
        }

    }

    private static boolean isSelfDescribingPlugin(String str) {

        try {

            String possibleClassName = str.substring(str.lastIndexOf('-') + 1);
            List<String> matches = Utils.getFullyQualifiedClassNames(possibleClassName);
            for (String match : matches) {
                try {
                    Class currentMatch = Class.forName(match);
                    Constructor constructor = currentMatch.getConstructor(Frame.class);
                    return true;
                } catch (Exception e) {
                    // do nothing
                }
            }

            try {
                Class possibleClass = Class.forName(possibleClassName);
                Constructor constructor = possibleClass.getConstructor(Frame.class);
                return true;
            } catch (Exception e) {
                // do nothing
            }

            return false;

        } catch (Exception e) {
            return false;
        }

    }

    private static Element createTag(Document doc, Element element, String tag) {
        String str = tag.substring(tag.lastIndexOf('-') + 1);
        Element tagElement = doc.createElement(str);
        element.appendChild(tagElement);
        return tagElement;
    }

    private static int createString(Document doc, Element element, String[] args, int index) throws IOException {
        String current = args[index].substring(args[index].lastIndexOf('-') + 1);
        if (args[index].startsWith("-")) {
            Element newElement = createTag(doc, element, current);
            while (true) {
                index++;
                if (index >= args.length) {
                    return index - 1;
                }
                if (args[index].startsWith("-")) {
                    return index - 1;
                }
                newElement.appendChild(doc.createTextNode(args[index]));
            }
        } else {
            element.appendChild(doc.createTextNode(current));
            return index;
        }
    }

    public static String[] readXMLAsArgs(String filename) {

        List temp = new ArrayList();

        try {

            File fXmlFile = new File(filename);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);
            doc.getDocumentElement().normalize();

            Element rootElement = doc.getDocumentElement();
            if (!(rootElement.getNodeName().equalsIgnoreCase("TasselPipeline"))) {
                throw new IllegalArgumentException("ReadPipelineXML: readXML: Root Node must be TasselPipeline: " + rootElement.getNodeName());
            }

            NodeList children = rootElement.getChildNodes();
            for (int i = 0; i < children.getLength(); i++) {
                Node current = children.item(i);
                getFlags(current, temp);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

        String[] args = new String[temp.size()];
        temp.toArray(args);
        return args;

    }

    private static void getFlags(Node node, List flags) {

        if (node.getNodeType() != Node.ELEMENT_NODE) {
            return;
        }

        String flagName = node.getNodeName().trim();
        flags.add("-" + flagName);
        NodeList children = node.getChildNodes();
        for (int i = 0; i < children.getLength(); i++) {
            Node current = children.item(i);
            if (current.getNodeType() == Node.TEXT_NODE) {
                String temp = current.getNodeValue().trim();
                if (temp.length() != 0) {
                    flags.add(temp);
                }
            } else {
                getFlags(current, flags);
            }
        }

        if (isSelfDescribingPlugin(flagName)) {
            flags.add("-endPlugin");
        }

    }
}
