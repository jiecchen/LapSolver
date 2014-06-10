/**
 * @file BindingGenerator.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Wed Jun 4 2014
 *
 * A Python code generator to produce bindings to packages within
 * jar files. Used to create the bindings for LapSolver.
 *
 * Depends on numpy and jpype, both of which are available
 * in Linux distribution package managers, and via pip.
 *
 */

package pybinding;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.*;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

public class BindingGenerator {
    private Map<String, Class> classNames;
    private List<String> packages;
    private final String topLevelPackage;
    private final String jarPath;

    /**
     * Construct a new BindingGenerator
     *
     * @param jarPath     path to the JAR for which to generate python bindings
     * @param packageName the name of the Java package (in the JAR) to be bound
     * @throws IOException when the JAR file cannot be opened.
     */
    public BindingGenerator(String jarPath, String packageName) throws IOException {
        this.topLevelPackage = packageName;
        this.jarPath = jarPath;

        File inputJar = new File(jarPath);
        JarFile module = new JarFile(inputJar);
        URLClassLoader jarLoader = URLClassLoader.newInstance(new URL[]{inputJar.toURI().toURL()});

        classNames = new HashMap<>();
        packages = new ArrayList<>();

        Enumeration<JarEntry> entries = module.entries();
        while (entries.hasMoreElements()) {
            JarEntry entry = entries.nextElement();
            String entryName = entry.getName().replace('/', '.');

            if (entryName.startsWith(packageName)) { // if the jar entry is in the package
                if (entryName.endsWith(".class")) { // if it's a class
                    String className = entryName.substring(0, entryName.length() - ".class".length());
                    if (className.contains("$")) { // TODO: support inner classes
                        System.err.println("Warning: skipping inner class `" + className + "`");
                    } else {
                        try { // attempt to load the class
                            classNames.put(className, jarLoader.loadClass(className));
                        } catch (ClassNotFoundException e) {
                            System.err.println("Warning: unable to load class `" + className + "`");
                        }
                    }
                } else if (entryName.endsWith(".")) { // if it's a package
                    packages.add(entryName.substring(0, entryName.length() - 1));
                } else {
                    System.err.println("Warning: unknown object `" + entryName + "`");
                }
            }
        }

        /*
         * Sorting by length ensures that a subpackage always comes after
         * its parent. This is true since the parent's name is a prefix
         * of the child's name.
         *
         * e.g. `pkgA.pkgB` will come after `pkgA`
         */
        Collections.sort(packages, new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                return o1.length() - o2.length();
            }
        });
    }

    /**
     * Find the containing package of a class name
     *
     * @param objectName fully-qualified class name (e.g. "java.lang.System")
     * @return Name of the package (e.g. "java.lang")
     */
    private String getContainingPackage(String objectName) {
        int lastIdx = objectName.lastIndexOf('.');
        if (lastIdx == -1)
            return "";
        return objectName.substring(0, lastIdx);
    }

    /**
     * Load a template string from resources
     *
     * @param templateFile the name of the template (e.g. "tlheader")
     * @return The contents of the template file
     */
    private String getTemplate(String templateFile) {
        InputStream headerInput = getClass().getResourceAsStream("resources/templates/" + templateFile + ".txt");
        Scanner reader = new Scanner(headerInput).useDelimiter("\\A");
        return reader.hasNext() ? reader.next() : "";
    }

    /**
     * Create a python-ready reference to a package or class
     *
     * @param packageName the name of the java package (e.g. "java.lang")
     * @return The jpype call to load it (e.g. "JPackage(\"java\").lang")
     */
    private String asJPype(String packageName) {
        String[] components = packageName.split("\\.");
        if (components.length == 0)
            return "";
        StringBuilder newName = new StringBuilder();
        newName.append("JPackage(\"");
        newName.append(components[0]);
        newName.append("\")");
        for (int i = 1; i < components.length; i++) {
            newName.append(".");
            newName.append(components[i]);
        }
        return newName.toString();
    }

    /**
     * Generate the python module.
     *
     * @param directoryName directory to store the generated module
     */
    public void generate(String directoryName) {
        String jarFileName = new File(jarPath).getName();

        // Make sure the directory name has consistent separators
        // and is clearly either absolute or relative
        String outputDirectory = directoryName;
        if (!outputDirectory.startsWith("/") && !outputDirectory.startsWith("./"))
            outputDirectory = "./" + outputDirectory;
        if (!outputDirectory.endsWith("/"))
            outputDirectory += "/";

        for (String curPackage : packages) {
            String dirName = outputDirectory + curPackage.replace('.', '/');
            File moduleFile = new File(dirName + "/__init__.py");

            // Create the package directory
            if (!moduleFile.getParentFile().mkdirs() && !moduleFile.getParentFile().exists()) {
                System.err.println("Error: could not create directory `" + dirName + "`");
                System.exit(0);
            }

            StringBuilder pythonModule = new StringBuilder();

            String header = getTemplate(curPackage.equals(topLevelPackage) ? "tlheader" : "subheader");
            header = header.replaceAll("\\{\\{JAR_FILE\\}\\}", jarFileName);
            header = header.replaceAll("\\{\\{PACKAGE\\}\\}", topLevelPackage);

            pythonModule.append(header).append("\n");

            for (Class currentClass : classNames.values()) {
                String className = currentClass.getName();
                String parentPackage = getContainingPackage(className);

                if (parentPackage.equals(curPackage)) {
                    if (currentClass.isInterface()) {
                        System.err.println("Warning: skipping interface " + className);
                        continue;
                    }

                    // The name Python will use to refer to the generated class
                    String pyClassName = className.substring(className.lastIndexOf('.') + 1);

                    // Emit the class declaration
                    pythonModule.append("class ").append(pyClassName).append("(_GeneratedObject):\n");

                    // Emit the constructor for the class (if it has one)
                    if (currentClass.getConstructors().length > 0) {
                        String ctor = getTemplate("constructor");
                        ctor = ctor.replaceAll("\\{\\{CLASS\\}\\}", pyClassName);
                        ctor = ctor.replaceAll("\\{\\{JAVA_PACKAGE\\}\\}", curPackage);
                        ctor = ctor.replaceAll("\\{\\{PYTHON_PACKAGE\\}\\}", asJPype(className));
                        pythonModule.append(ctor).append("\n");
                    } else {
                        System.err.println("Warning: no constructors found for " + className);
                    }

                    // Generate wrappers for static methods.
                    for (Method method : currentClass.getMethods()) {
                        if(Modifier.isStatic(method.getModifiers())) {
                            String met = getTemplate("staticmethod");
                            met = met.replaceAll("\\{\\{CLASS\\}\\}", pyClassName);
                            met = met.replaceAll("\\{\\{METHOD\\}\\}", method.getName());
                            pythonModule.append(met).append("\n");
                        }
                    }

                    // TODO: generate code to properly handle fields

                    pythonModule.append("\n\n");
                }
            }

            // Write the generated module to __init__.py.
            try (BufferedWriter writer = Files.newBufferedWriter(moduleFile.toPath(),
                    Charset.forName("US-ASCII"))) {
                String contents = pythonModule.toString();
                writer.write(contents, 0, contents.length());
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(0);
            }
        }

        // Copy the jar file to the new directory
        try {
            File destination = new File(outputDirectory + topLevelPackage + "/_src/" + jarFileName);
            destination.mkdirs();
            Files.copy(new File(jarPath).toPath(), destination.toPath(), REPLACE_EXISTING);
        } catch (IOException e) {
            System.err.println("Error: could not copy JAR file to Python module directory.");
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        try {
            BindingGenerator bindingGenerator = new BindingGenerator(args[0], args[1]);
            bindingGenerator.generate(args[2]);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
