package lapsolver;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;

public class Logger {

    public boolean opened;
    FileWriter fw;
    String name;

    public Logger() {
        opened = false;
    }

    public void start(String name) {

        this.name = name;

        try {
            fw = new FileWriter(name);
            Date date = new Date();
            fw.write("Logger " + name + " at " + date.toString() + "\n");
            fw.flush();
        } catch (IOException e) {
            throw new Error("Error writing to " + name);
        }

        opened = true;
    }

    public void write(String str) {
        if (opened) {
            try {
                fw.write(str + "\n");
                fw.flush();
            } catch (IOException e) {
                throw new Error("Error writing to " + name);
            }
        }
    }

    public void flush() {
        if (opened) {
            try {
                fw.flush();
            } catch (IOException e) {
                throw new Error("Error writing to " + name);
            }
        }
    }

    public void done() {
        if (opened) {
            try {
                fw.close();
            } catch (IOException e) {
                throw new Error("Error writing to " + name);
            }
        }
        opened = false;
    }
}
