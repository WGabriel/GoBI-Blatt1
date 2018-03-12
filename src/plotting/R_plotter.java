package plotting;

import java.io.*;

public class R_plotter implements Runnable {

    public static void main(String[] args) {
        System.out.println("Main started.");
        InputStream stdout = null;
        InputStream stderr = null;

        try {
            String rCom = "print(\'Hi\');";
            Process p = new ProcessBuilder("R", "-e", rCom).start();
            p.waitFor();
            stdout = p.getInputStream();
            stderr = p.getErrorStream();
            BufferedReader bri = new BufferedReader(new InputStreamReader(stdout));

            String line;
            while ((line = bri.readLine()) != null) {
                System.out.println("Line: " + line);
            }


        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public void run() {
    }
}
