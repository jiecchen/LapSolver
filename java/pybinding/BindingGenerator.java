package pybinding;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.*;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

/**
 * Created by alex on 6/4/14.
 */
public class BindingGenerator {
    private Map<String, Class> classNames;
    private List<String> packages;
    private String topLevelPackage;
    private String jarPath;

    public BindingGenerator(String jarPath, String packageName) throws IOException {
        this.topLevelPackage = packageName;
        this.jarPath = jarPath;

        File file = new File(jarPath);
        JarFile archive = new JarFile(file);
        URLClassLoader loader = URLClassLoader.newInstance(new URL[]{file.toURI().toURL()});

        classNames = new HashMap<>();
        packages = new ArrayList<>();

        Enumeration<JarEntry> entries = archive.entries();
        while(entries.hasMoreElements()) {
            JarEntry entry = entries.nextElement();
            String objectName = entry.getName().replace('/', '.');

            if (objectName.startsWith(packageName)) {
                if (objectName.endsWith(".class")) {
                    String className = objectName.substring(0, objectName.length() - ".class".length());
                    if (className.contains("$")) {
                        System.err.println("Warning: skipping inner class `" + className + "`");
                    } else {
                        try {
                            classNames.put(className, loader.loadClass(className));
                        } catch (ClassNotFoundException e) {
                            System.err.println("Warning: unable to load class `" + className + "`");
                        }
                    }
                } else if (objectName.endsWith(".")) {
                    packages.add(objectName.substring(0, objectName.length() - 1));
                } else {
                    System.err.println("Warning: unknown object `" + objectName + "`");
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

    private String getContainingPackage(String objectName) {
        int lastIdx = objectName.lastIndexOf('.');
        if(lastIdx == -1)
            return "";
        return objectName.substring(0,lastIdx);
    }

    private String getTemplate(String templateFile) {
        InputStream headerInput = getClass().getResourceAsStream("resources/templates/" + templateFile + ".txt");
        Scanner reader = new Scanner(headerInput).useDelimiter("\\A");
        return reader.hasNext() ? reader.next() : "";
    }

    private String asJPype(String packageName) {
        String[] components = packageName.split("\\.");
        if(components.length == 0)
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

    public void generate(String outputDirectory) {
        String jarFileName = new File(jarPath).getName();

        String prefix = outputDirectory;
        if(!prefix.startsWith("/") && !prefix.startsWith("./"))
            prefix = "./" + prefix;
        if(!prefix.endsWith("/"))
            prefix += "/";

        // Create the package directory
        for (String curPackage : packages) {
            String dirName = prefix + curPackage.replace('.', '/');
            File moduleFile = new File(dirName + "/__init__.py");
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
                if(parentPackage.equals(curPackage)) {
                    if(currentClass.isInterface()) {
                        System.err.println("Warning: skipping interface " + className);
                        continue;
                    }

                    String pyClassName = className.substring(className.lastIndexOf('.') + 1);
                    pythonModule.append("class ").append(pyClassName).append("(_GeneratedObject):\n");
                    if(currentClass.getConstructors().length > 0) {
                        String ctor = getTemplate("constructor");
                        ctor = ctor.replaceAll("\\{\\{CLASS\\}\\}", pyClassName);
                        ctor = ctor.replaceAll("\\{\\{JAVA_PACKAGE\\}\\}", curPackage);
                        ctor = ctor.replaceAll("\\{\\{PYTHON_PACKAGE\\}\\}", asJPype(className));
                        pythonModule.append(ctor).append("\n");
                    } else {
                        System.err.println("Warning: no constructors found for " + className);
                    }
                    pythonModule.append("\n\n"); // obey PEP8
                }
            }

            try(BufferedWriter writer = Files.newBufferedWriter(moduleFile.toPath(),
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
            File destination = new File(prefix + topLevelPackage + "/_src/" + jarFileName);
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
            bindingGenerator.generate("python");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
