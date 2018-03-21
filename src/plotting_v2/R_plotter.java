package plotting_v2;

import java.io.*;
import java.util.*;

public class R_plotter implements Runnable {

    public static void main(String[] args) {
        // ATTENTION: To swap between ExonSkipping<->BaseSkipping: modify the values "isExon", "xlim" (R) and "xaxp" (R)
        System.out.println("Main started.");

        //Add all files in output to ArrayList "outputs"
        File path = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt1\\plots-v2\\outputs");
        ArrayList<File> outputs = new ArrayList<>(Arrays.asList(Objects.requireNonNull(path.listFiles())));

        // ---Gather information for the plots---
        // for exonSkipping -> isExon = true. for baseSkipping -> isExon = false.
        boolean isExon = true;
        TreeMap<Integer, Event> allSkippings = new TreeMap<>(Collections.reverseOrder());
        for (File out : outputs) {
            // Key = #skippedExons/skippedBases, value = Event
            TreeMap<Integer, Event> maxSkippings = getSkippingData(out, isExon);
            for (Map.Entry<Integer, Event> entry : maxSkippings.entrySet()) {
                if (allSkippings.containsKey(entry.getKey())) {
                    allSkippings.get(entry.getKey()).geneIds.addAll(entry.getValue().geneIds);
                    allSkippings.get(entry.getKey()).numberOfEvents += entry.getValue().numberOfEvents;
                } else
                    allSkippings.put(entry.getKey(), entry.getValue());
            }
        }
        System.out.println("The Top genes with most skippings are listed below. isExon: " + isExon);
        for (Map.Entry<Integer, Event> entry : allSkippings.entrySet()) {
            //System.out.println("Skippings: " + entry.getKey() + " Event: " + entry.getValue().numberOfEvents + " " + entry.getValue().geneIds.toString());
        }
        // ---END: Gather information for the plots---

        // ###########RESULTS###########
        // Max Events: 17424
        // Max ExonSkipping: 169
        // Max BaseSkipping: 26106
        // ###########RESULTS###########


        // ---Create and execute R-Command---
        String rCom = "png(\"C:/Users/Gabriel/Desktop/GoBI/Blatt1/plots-v2/R/test.png\", width = 1200, height = 548)\n"
                + "plot(1, 1, type = \"n\", xlim=c(0,170),xaxp=c(0,170,17),ylim=c(0,17424),yaxp=c(0,17500,7),ann=F,panel.first=grid())\n"
                + "title(main=\"Empirical Cumulative Distribution Function\",xlab=\"Max. number of Skipped Bases\",ylab=\"Cumulative number of events\")\n";

        ArrayList<String> colors = new ArrayList<>(Arrays.asList("pink", "blue", "orange", "violet", "green", "red"));
        int colorcounter = 0;
        for (File out : outputs) {
            // Key = #skippedExons/skippedBases, value = Event
            TreeMap<Integer, Event> maxSkippingsReverse = new TreeMap<>();
            maxSkippingsReverse.putAll(getSkippingData(out, isExon));

            int numberOfEvents = 0;
            int maxSkipped = 0;
            for (Map.Entry<Integer, Event> entry : maxSkippingsReverse.entrySet()) {
                maxSkipped = entry.getKey();
                numberOfEvents += entry.getValue().numberOfEvents;
                rCom += "points(" + maxSkipped + "," + numberOfEvents + ", col = \"" + colors.get(colorcounter) + "\", pch=\".\")\n";
                //System.out.println("Cumul. Skippings: " + maxSkipped + " Cumul. NrOfEvents: " + numberOfEvents + " Genes: " + entry.getValue().geneIds.toString());
            }
            colorcounter++;
        }

        rCom += "legend(\"bottomright\","
                + "title=\"GTF source\","
                + "legend=c(\"gencode.v25\",\"Homo_sapiens.GRCh37.67\",\"Homo_sapiens.GRCh37.75\",\"Homo_sapiens.GRCh38.86\",\"Homo_sapiens.GRCh38.90\",\"Mus_musculus.GRCm38.75\"),"
                + "pch=c(16,16,16,16,16,16),"
                + "col=c(\"pink\", \"blue\",\"orange\",\"violet\",\"green\",\"red\"))\n"
                + "dev.off()\n";
        //String rCom = "print(\'Hi\');\n";
        File rScript = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt1\\plots-v2\\R\\myScript.R");
        writeStringToFile(rCom, rScript);
        runRscript(rScript);
        // ---END: Create and execute R-Command---
    }

    private static void runRscript(File script) {
        try {
            System.out.println("R-Script: " + script.getAbsolutePath());
            Process p = new ProcessBuilder("Rscript", "\"" + script.getAbsolutePath() + "\"").start();
            p.waitFor();
            String line;
            // --- stdout ---
            BufferedReader stdout = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((line = stdout.readLine()) != null)
                System.out.println("StdOut: " + line);
            // --- sterr ---
            BufferedReader stderr = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            while ((line = stderr.readLine()) != null)
                System.err.println("StdErr: " + line);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void writeStringToFile(String string, File file) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(file));
            writer.write(string);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static TreeMap<Integer, Event> getSkippingData(File file, boolean isExon) {
        // return either exonSkippings (->exon=true) or baseSkippings (->exon=false)
        System.out.println("Begin: Parsing SkippingData. isExon: " + isExon + ". File: " + file.getAbsolutePath());
        // key = #maxExonSkipping / #maxBaseSkipping, value = gene_id's
        int lineNumber = 0; // for maximum plot ranges
        TreeMap<Integer, Event> result = new TreeMap<>(Collections.reverseOrder());
        int column;
        if (isExon)
            column = 11; //max_exons in column 11
        else
            column = 13; //max_bases in columt 13
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = br.readLine(); //skip Headers
            while ((line = br.readLine()) != null) {
                lineNumber++;
                if (line.split("\t").length < 14) {
                    System.err.println("14 tabs required! Found: " + line); // error checking
                } else {
                    String gene_id = line.split("\t")[0];
                    Integer max_skipped = Integer.parseInt(line.split("\t")[column]);

                    if (result.containsKey(max_skipped)) {
                        result.get(max_skipped).geneIds.add(gene_id);
                        result.get(max_skipped).numberOfEvents++;
                    } else {
                        Event event = new Event(new HashSet<>(Collections.singletonList(gene_id)), 1);
                        result.put(max_skipped, event);
                    }
                }
            }
            // System.out.println("LineNumber w/o header: " + lineNumber);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return result;
    }

    @Override
    public void run() {
    }

}

