package lapsolver.util;

import java.io.*;
import java.util.HashSet;

public class NativeLoader {
    private NativeLoader() { }
    private static HashSet<String> loadedLibraries = new HashSet<>();

    public static void loadLibrary(String libName)
    {
        if(!loadedLibraries.contains(libName)) {
            try {
                System.loadLibrary(libName);
            } catch (UnsatisfiedLinkError e) {
                try {
                    String nativeName = System.mapLibraryName(libName);
                    int start = nativeName.indexOf(libName);
                    String prefix = nativeName.substring(0, start) + "jni";
                    String suffix = nativeName.substring(start + libName.length());

                    File tempLib = File.createTempFile(prefix, suffix);
                    tempLib.deleteOnExit();

                    if(!tempLib.exists())
                        throw new FileNotFoundException("File " + tempLib.getAbsolutePath() + " could not be created.");

                    InputStream in = NativeLoader.class.getResourceAsStream("/lib/" + nativeName);
                    if(in == null)
                        throw new FileNotFoundException("File /lib/" + nativeName + " was not found in the JAR");

                    try (OutputStream out = new FileOutputStream(tempLib)) {
                        byte[] buffer = new byte[1024];
                        int readBytes;

                        while ((readBytes = in.read(buffer)) != -1) {
                            out.write(buffer, 0, readBytes);
                        }
                    } finally {
                        in.close();
                    }

                    // Finally, load the library
                    System.load(tempLib.getAbsolutePath());
                } catch (IOException e1) {
                    throw new RuntimeException(e1);
                }
            }
            loadedLibraries.add(libName);
        }
    }
}
