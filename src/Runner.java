import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class Runner {

    public static void main(String[] args) {
        // Check if console parameters are correct
        if (args.length != 4) {
            System.err.println("Too few or many ConsoleParameters.");
            System.exit(1);
        }
        File gtfInput = null;
        File tsvOutput = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-o")) {
                tsvOutput = new File(args[i + 1]);
            }
            if (args[i].equals("-gtf")) {
                gtfInput = new File(args[i + 1]);
            }
        }
        assert gtfInput != null;
        assert tsvOutput != null;
        System.out.println("gtf input file:\t" + gtfInput.getAbsolutePath());
        System.out.println("tsv output file:\t" + tsvOutput.getAbsolutePath());

        HashMap<String, Gene> allGenes = parseGtf(gtfInput);
        findIntrons(allGenes);

        HashSet<OutputLine> outputLines = getOutput(allGenes);

        writeOutput(outputLines, tsvOutput);

        System.out.println("End of main-method.");
    }

    private static void findIntrons(HashMap<String, Gene> allGenes) {
        System.out.println("|findIntrons| Started...");
        for (Gene gene : allGenes.values()) {
            for (Protein protein : gene.proteins.values()) {
                //Go through all proteins
                CDS previousCds = null;
                for (CDS currentCDS : protein.cdss) {
                    // Go through all CDS
                    if (previousCds != null) {
                        if (currentCDS.start < previousCds.end) {
                            protein.introns.add(new Intron(previousCds.end, currentCDS.start));
                            System.out.println("intron added");
                        } else {
                            System.err.println("|findIntrons| Intron start > Intron end. In CDS: " + currentCDS.start + "-" + currentCDS.end + " (" + currentCDS.protein_id + ")");
                        }
                        previousCds = new CDS(currentCDS);
                    }
                }
            }
        }
        System.out.println("|findIntrons| Finished.");
    }

    private static HashSet<OutputLine> getOutput(HashMap<String, Gene> allGenes) {
        HashSet<OutputLine> result = new HashSet<>();
        for (Gene gene : allGenes.values()) {
            for (Protein protein : gene.proteins.values()) {
                for (Intron intr : protein.introns) {
                    // For every Intron
                    OutputLine l = new OutputLine();
                    CDS anyCDS = protein.cdss.iterator().next();
                    l.gene_id = anyCDS.gene_id;
                    l.gene_symbol = anyCDS.gene_name;
                    l.chromosome = anyCDS.chr;
                    l.strand = anyCDS.strand;
                    l.nprots = gene.proteins.size();
                    l.ntrans = gene.transcripts.size();
                    l.sv = new Intron(intr.start, intr.end);
                    l.wt = new TreeSet<>(new Intron()); // WT introns within the SV intron separated by | as start:end
                    l.sv_prots = new HashSet<>(); // ids of the SV CDS-s, separated by |
                    l.wt_prots = new HashSet<>(); // ids of the WT CDS-s, separated by |
                    l.min_skipped_exon = Integer.MAX_VALUE;
                    l.max_skipped_exon = Integer.MIN_VALUE;
                    l.min_skipped_base = Integer.MAX_VALUE;
                    l.max_skipped_base = Integer.MIN_VALUE;

                    for (Protein tempProtein : gene.proteins.values()) {
                        // find leftIntronBorder, rightIntronBorder, withinIntronRange
                        HashSet<CDS> leftIntronBorder = new HashSet<>();
                        HashSet<CDS> rightIntronBorder = new HashSet<>();
                        HashSet<CDS> withinIntronRange = new HashSet<>();
                        if (!tempProtein.protein_id.equals(protein.protein_id)) {
                            for (CDS cds : tempProtein.cdss) {
                                if (cds.end == intr.start)
                                    leftIntronBorder.add(new CDS(cds));
                                else if (cds.start == intr.end)
                                    rightIntronBorder.add(new CDS(cds));
                                else if (cds.start >= intr.start && cds.end <= intr.end)
                                    withinIntronRange.add(new CDS(cds));
                            }
                        }

                        //skipped CDS -> Protein-Intersection in lefIntronBorder, rightIntronBorder and withinIntronRange
                        HashSet<CDS> skippedCDS = new HashSet<>(leftIntronBorder);
                        skippedCDS.retainAll(rightIntronBorder);
                        skippedCDS.retainAll(withinIntronRange);

                        if (!skippedCDS.isEmpty()) {
                            if (skippedCDS.size() < l.min_skipped_exon)
                                l.min_skipped_exon = skippedCDS.size();
                            if (skippedCDS.size() > l.max_skipped_exon)
                                l.max_skipped_exon = skippedCDS.size();

                            // joint size of all skipped CDS
                            int jointLengthCDS = 0;
                            CDS previousSkippedCDS = null;
                            for (CDS aSkippedCDS : skippedCDS) {
                                jointLengthCDS += (aSkippedCDS.end - aSkippedCDS.start);
                                l.wt_prots.add(aSkippedCDS.protein_id);

                                if (previousSkippedCDS != null) {
                                    l.wt.add(new Intron(previousSkippedCDS.end, aSkippedCDS.start));
                                }
                                previousSkippedCDS = aSkippedCDS;
                            }

                            if (jointLengthCDS < l.min_skipped_base)
                                l.min_skipped_base = jointLengthCDS;
                            if (jointLengthCDS > l.max_skipped_base)
                                l.max_skipped_base = jointLengthCDS;

                            // gather sv_prots
                            HashSet<CDS> sv_prots = new HashSet<>(leftIntronBorder);
                            sv_prots.retainAll(rightIntronBorder);
                            sv_prots.removeAll(skippedCDS);
                            for (CDS cds : sv_prots) {
                                l.sv_prots.add(cds.protein_id);
                            }

                            // error checking
                            if (l.gene_id.isEmpty() || l.gene_symbol.isEmpty() || l.chromosome.isEmpty() || l.strand.isEmpty()
                                    || l.ntrans == 0 || l.ntrans == 0 || l.sv == null || l.wt.isEmpty() || l.sv_prots.isEmpty()
                                    || l.wt_prots.isEmpty() || l.min_skipped_base < 1 || l.max_skipped_base < 1 || l.min_skipped_exon < 1
                                    || l.max_skipped_exon < 1) {
                                System.err.println("|getOutput| One value in OutputLine is empty. Intron: " + intr.start + "-" + intr.end);
                            }
                            result.add(l);
                        }
                    }
                }
            }
        }
        return result;
    }


    public static HashMap<String, Gene> parseGtf(File gtfInput) {
        // Key is the always the respective ID
        HashMap<String, Gene> allGenes = new HashMap<>();
        int geneCounter = 0;
        int proteinCounter = 0;
        int cdsCounter = 0;
        // read gtfInput-file
        try {
            BufferedReader br = new BufferedReader(new FileReader(gtfInput));
            String line;
            int linecounter = 0; // for error printing only
            System.out.println("Begin: Parsing file.");
            while ((line = br.readLine()) != null) {
                linecounter++;
                // ignore comments (beginning with "#")
                if (!line.startsWith("#")) {
                    // For every line in input
                    String[] tabSeparated = line.split("\\t");
                    String seqname = tabSeparated[0];
                    // String source = tabSeparated[1];
                    String feature = tabSeparated[2];
                    int start = Integer.parseInt(tabSeparated[3]);
                    int end = Integer.parseInt(tabSeparated[4]);
                    // String score = tabSeparated[5];
                    String strand = tabSeparated[6];
                    // String frame = tabSeparated[7];

                    String attribute = tabSeparated[8];
                    // parameters needed to construct new CDS
                    String protein_id = "";
                    String transcript_id = "";
                    String transcript_name = "";
                    String gene_id = "";
                    String gene_name = "";
                    // gather parameters from String "attribute"
                    String[] attributeSeparated = attribute.split(";");
                    // search in attributeSeparated for parameters
                    for (String attr : attributeSeparated) {
                        if (attr.contains("protein_id")) {
                            // get only value between quotation marks
                            protein_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("transcript_id")) {
                            transcript_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("transcript_name")) {
                            transcript_name = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("gene_id")) {
                            gene_id = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        } else if (attr.contains("gene_name")) {
                            gene_name = attr.substring(attr.indexOf("\"") + 1, attr.lastIndexOf("\""));
                        }
                    }

                    // -------For lines, which are CDS:-------
                    if (feature.equalsIgnoreCase("CDS")) {
                        // Construct cds and add to exonList
                        CDS cds = new CDS(seqname, start, end, strand, protein_id, transcript_id, transcript_name, gene_id, gene_name);
                        // Check, if all cds values are actually filled
                        if (cds.start == 0 || cds.end == 0 || cds.chr.isEmpty() || cds.strand.isEmpty() || cds.protein_id.isEmpty() || cds.transcript_id.isEmpty() || cds.gene_id.isEmpty()) {
                            System.err.println("CDS in line " + linecounter + " has an empty value!");
                        }
                        // Create data structure allGenes, contains proteins, contains cdss
                        if (!allGenes.containsKey(cds.gene_id)) {
                            // add a new gene to allGenes, which first requires a set of Proteins, which requires a set of CDS
                            TreeSet<CDS> cdss = new TreeSet<>(new CDS());
                            cdss.add(new CDS(cds));

                            HashMap<String, Protein> proteins = new HashMap<>();
                            proteins.put(cds.protein_id, new Protein(cds.protein_id, cdss));

                            allGenes.put(cds.gene_id, new Gene(cds.gene_id, proteins));
                            geneCounter++;
                            proteinCounter++;
                            cdsCounter++;

                            allGenes.get(cds.gene_id).transcripts.add(cds.transcript_id);
                        } else {
                            // check, if existing gene already contains the found transcript
                            Gene currentGene = allGenes.get(cds.gene_id);
                            if (!currentGene.proteins.containsKey(cds.protein_id)) {
                                // add a new protein to proteins, which first requires a set of cdss
                                TreeSet<CDS> cdss = new TreeSet<>(new CDS());
                                cdss.add(new CDS(cds));

                                currentGene.proteins.put(cds.protein_id, new Protein(cds.protein_id, cdss));
                                proteinCounter++;
                                cdsCounter++;

                                allGenes.get(cds.gene_id).transcripts.add(cds.transcript_id);
                            } else {
                                // If protein already exists, just add new CDS to this proteins
                                currentGene.proteins.get(cds.protein_id).cdss.add(new CDS(cds));
                                cdsCounter++;

                                allGenes.get(cds.gene_id).transcripts.add(cds.transcript_id);
                            }
                        }
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Finished: Parsing file.");
        System.out.println("\tcdsCounter: " + cdsCounter);
        System.out.println("\tprotein:" + proteinCounter);
        System.out.println("\tgeneCounter:" + geneCounter);
        return allGenes;
    }

    private static void writeOutput(HashSet<OutputLine> outputLines, File tsvOutput) {
        System.out.println("|writeOutput| Begin...");
        try {
            PrintWriter out = new PrintWriter(tsvOutput, "UTF-8");
            out.println("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases");
            for (OutputLine l : outputLines) {
                // Convert the HashSet<String> wt_prots to String
                String wt_prots_string = "";
                if (l.wt_prots != null) {
                    for (String s : l.wt_prots) {
                        wt_prots_string = wt_prots_string.concat(s).concat("|");
                    }
                    // crop last "|"
                    if (wt_prots_string.endsWith("|")) {
                        wt_prots_string = wt_prots_string.substring(0, wt_prots_string.length() - 1);
                    }
                }

                // Convert the HashSet<String> sv_prots to String
                String sv_prots_string = "";
                if (l.sv_prots != null) {
                    for (String s : l.sv_prots) {
                        sv_prots_string = sv_prots_string.concat(s).concat("|");
                    }
                    // crop last "|"
                    if (sv_prots_string.endsWith("|")) {
                        sv_prots_string = sv_prots_string.substring(0, sv_prots_string.length() - 1);
                    }
                }

                // Convert the HashSet<String> wt_prots to String
                String wt_string = "";
                if (l.wt != null) {
                    wt_string = Integer.toString(l.sv.start);
                    for (Intron intr : l.wt) {
                        wt_prots_string += ":" + Integer.toString(intr.start) + "|" + Integer.toString(intr.end);
                    }
                    wt_string += ":" + Integer.toString(l.sv.end);
                }

                String thisLine = l.gene_id + "\t" + l.gene_symbol + "\t" + l.chromosome + "\t" + l.strand
                        + "\t" + l.nprots + "\t" + l.ntrans + "\t" + (l.sv.start + 1) + ":" + l.sv.end + "\t"
                        + wt_string + "\t" + wt_prots_string + "\t" + sv_prots_string + "\t"
                        + l.min_skipped_exon + "\t" + l.max_skipped_exon + "\t" + (l.min_skipped_base + 1)
                        + "\t" + (l.max_skipped_base + 1);
                out.println(thisLine);
                //System.out.println(thisLine);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("|writeOutput| Finished.");
    }

}
