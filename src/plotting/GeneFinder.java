package plotting;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

public class GeneFinder {
    // finds the top 10 genes for max ExonSkipping and BasePairs
    public static void main(String[] args) {
        File exonSkippingOutput = new File(
                "C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt1\\output\\Mus_musculus.GRCm38.75.tsv");

        findGenes(exonSkippingOutput);
    }

    private static void findGenes(File exonSkippingFile) {
        System.out.println("Begin: Parsing file.");
        // key = Gene_id, value = #exonSkippings / #baseSkippings
        HashMap<String, Integer> exonsSkippings = new HashMap<>();
        HashMap<String, Integer> baseSkippings = new HashMap<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(exonSkippingFile));
            String line = br.readLine(); //skip Headers
            while ((line = br.readLine()) != null) {
                if (line.split("\t").length < 14) {
                    // error checking
                    System.err.println("14 tabs required! Found: " + line);
                } else {
                    String gene_id = line.split("\t")[0];
                    Integer max_skipped_bases = Integer.parseInt(line.split("\t")[13]);
                    Integer max_skipped_exons = Integer.parseInt(line.split("\t")[11]);

                    if (exonsSkippings.containsKey(gene_id))
                        exonsSkippings.put(gene_id, exonsSkippings.get(gene_id) + max_skipped_exons);
                    else
                        exonsSkippings.put(gene_id, max_skipped_exons);

                    if (baseSkippings.containsKey(gene_id))
                        baseSkippings.put(gene_id, baseSkippings.get(gene_id) + max_skipped_bases);
                    else
                        baseSkippings.put(gene_id, max_skipped_bases);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        Map<String, Integer> sortedExons = sortByValues(exonsSkippings);
        Map<String, Integer> sortedBases = sortByValues(baseSkippings);

        System.out.println("---Most Exons:---");
        int counter = 1;
        for (Map.Entry<String, Integer> entry : sortedExons.entrySet()) {
            if (counter < 11) {
                // System.out.println(counter + ". Gene: " + entry.getKey() + " | #Exons: " + entry.getValue());
                String gene_id = entry.getKey();
                if (gene_id.contains(".")) {
                    // chops ".25" from ENSG00000155657.25
                    gene_id = gene_id.substring(0, gene_id.lastIndexOf('.'));
                }
                System.out.println("<li><a href=\"http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=" + gene_id + "\">" + gene_id + "(" + entry.getValue() + " Exons)</a></li>");
                counter++;
            }
        }
        System.out.println("---Most Bases:---");
        counter = 1;
        for (Map.Entry<String, Integer> entry : sortedBases.entrySet()) {
            if (counter < 11) {
                // System.out.println(counter + ": Gene: " + entry.getKey() + " | #Bases: " + entry.getValue());
                String gene_id = entry.getKey();
                if (gene_id.contains(".")) {
                    // chops ".25" from ENSG00000155657.25
                    gene_id = gene_id.substring(0, gene_id.lastIndexOf('.'));
                }
                System.out.println("<li><a href=\"http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=" + gene_id + "\">" + gene_id + "(" + entry.getValue() + " Bases)</a></li>");
                counter++;
            }
        }
        System.out.println("Finished: Parsing file.");
    }


    private static HashMap sortByValues(HashMap map) {
        // sorts a given hashmap by value (ascending, if value = Integer)
        List list = new LinkedList(map.entrySet());
        // Defined Custom Comparator here
        Collections.sort(list, new Comparator() {
            public int compare(Object o1, Object o2) {
                return ((Comparable) ((Map.Entry) (o2)).getValue()).compareTo(((Map.Entry) (o1)).getValue());
            }
        });

        // Here I am copying the sorted list in HashMap using LinkedHashMap to preserve the insertion order
        HashMap sortedHashMap = new LinkedHashMap();
        for (Iterator it = list.iterator(); it.hasNext(); ) {
            Map.Entry entry = (Map.Entry) it.next();
            sortedHashMap.put(entry.getKey(), entry.getValue());
        }
        return sortedHashMap;
    }

}
