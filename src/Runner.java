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
        int intronCounter = 0;
        for (Gene gene : allGenes.values()) {
            //System.out.println("Current gene: " + gene.gene_id);
            for (Protein protein : gene.proteins.values()) {
                // Go through all proteins
                // System.out.println("Current protein: " + protein.protein_id);
                CDS previousCds = null;
                for (CDS currentCDS : protein.cdss) {
                    // System.out.println("Current CDS: " + currentCDS.start + "-" + currentCDS.end);
                    // Go through all CDS
                    if (previousCds != null) {
                        if (previousCds.end < currentCDS.start) {
                            gene.introns.add(new Intron(previousCds.end, currentCDS.start));
                            //System.out.println("Intron added to Protein: " + protein.protein_id + ". Intr: " + previousCds.end + "-" + currentCDS.start);
                            intronCounter++;
                        } else {
                            System.err.println("|findIntrons| Intron start > Intron end. In CDS: " + previousCds.end + "-" + currentCDS.start + " (" + currentCDS.protein_id + ")");
                        }
                    }
                    previousCds = new CDS(currentCDS);
                }
            }
        }
        System.out.println("|findIntrons| Finished. Added " + intronCounter + " Introns.");
    }

    private static HashSet<OutputLine> getOutput(HashMap<String, Gene> allGenes) {
        System.out.println("|getOutput| Started...");
        HashSet<OutputLine> result = new HashSet<>();
        for (Gene gene : allGenes.values()) {
            // System.out.println("Going through introns in gene: " + gene.gene_id);
            for (Intron intr : gene.introns) {
                // For every Intron
                OutputLine l = new OutputLine();
                CDS anyCDS = gene.proteins.values().iterator().next().cdss.iterator().next();
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

                for (Protein protein : gene.proteins.values()) {
                    // System.out.println("\t...checking against tempProtein: " + tempProtein.protein_id + " containing " + tempProtein.cdss.size() + " CDS's.");
                    // find leftIntronBorder, rightIntronBorder, withinIntronRange
                    HashSet<CDS> leftIntronBorder = new HashSet<>();
                    HashSet<CDS> rightIntronBorder = new HashSet<>();
                    // contains skipped Exons
                    HashSet<CDS> withinIntronRange = new HashSet<>();
                    for (CDS cds : protein.cdss) {
                        if (cds.end == intr.start)
                            leftIntronBorder.add(new CDS(cds));
                        else if (cds.start == intr.end)
                            rightIntronBorder.add(new CDS(cds));
                        else if (cds.start >= intr.start && cds.end <= intr.end)
                            withinIntronRange.add(new CDS(cds));
                    }

                    if (!withinIntronRange.isEmpty() && !leftIntronBorder.isEmpty() && !rightIntronBorder.isEmpty()) {
                        // withinIntronRange--> skippedExons
                        if (withinIntronRange.size() < l.min_skipped_exon)
                            l.min_skipped_exon = withinIntronRange.size();
                        if (withinIntronRange.size() > l.max_skipped_exon)
                            l.max_skipped_exon = withinIntronRange.size();

                        // joint size of all skipped CDS
                        int jointLengthCDS = 0;
                        for (CDS aSkippedCDS : withinIntronRange) {
                            jointLengthCDS += (aSkippedCDS.end - aSkippedCDS.start);
                            l.wt_prots.add(aSkippedCDS.protein_id);
                            l.wt.add(new Intron(aSkippedCDS.start, aSkippedCDS.end));
                        }

                        if (jointLengthCDS < l.min_skipped_base)
                            l.min_skipped_base = jointLengthCDS;
                        if (jointLengthCDS > l.max_skipped_base)
                            l.max_skipped_base = jointLengthCDS;

                    } else if (withinIntronRange.isEmpty() && !leftIntronBorder.isEmpty() && !rightIntronBorder.isEmpty()) {
                        // gather sv_prots
                        for (CDS cds : leftIntronBorder)
                            l.sv_prots.add(cds.protein_id);
                        for (CDS cds : rightIntronBorder)
                            l.sv_prots.add(cds.protein_id);
                    }
                }
                if (l.max_skipped_exon > 0) {
                    //
                    //
                    // System.out.println("added line: "+getOutputString(l));
                    // error checking
                    checkIfOutputLineIsEmpty(l);
                    result.add(l);
                }
            }
        }
        System.out.println("|getOutput| Finished. Added " + result.size() + " result lines.");
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

                    // Update ntrans
                    if (!transcript_id.isEmpty() && !gene_id.isEmpty()) {
                        if (!allGenes.containsKey(gene_id)) {
                            allGenes.put(gene_id, new Gene(gene_id));
                            geneCounter++;
                        }
                        allGenes.get(gene_id).transcripts.add(transcript_id);
                    }

                    // -------For lines, which are CDS:-------
                    if (feature.equalsIgnoreCase("CDS")) {
                        // Construct cds and add to exonList
                        CDS cds = new CDS(seqname, start, end, strand, protein_id, transcript_id, transcript_name, gene_id, gene_name);
                        // Check, if all cds values are actually filled
                        if (cds.chr.isEmpty() || cds.start < 1 || cds.end < 1 || cds.strand.isEmpty() || cds.protein_id.isEmpty() || cds.transcript_id.isEmpty() || cds.gene_id.isEmpty()) {
                            System.err.println("CDS in line " + linecounter + " has an empty value!");
                        }
                        // Create data structure allGenes, contains proteins, contains cdss
                        // check, if existing gene already contains the found transcript
                        Gene currentGene = allGenes.get(cds.gene_id);
                        if (!currentGene.proteins.containsKey(cds.protein_id)) {
                            // add a new protein to proteins, which first requires a set of cdss
                            TreeSet<CDS> cdss = new TreeSet<>(new CDS());
                            cdss.add(new CDS(cds));

                            currentGene.proteins.put(cds.protein_id, new Protein(cds.protein_id, cdss));
                            proteinCounter++;
                            cdsCounter++;
                        } else {
                            // If protein already exists, just add new CDS to this proteins
                            currentGene.proteins.get(cds.protein_id).cdss.add(new CDS(cds));
                            cdsCounter++;
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
                out.println(getOutputString(l));
                //System.out.println(thisLine);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("|writeOutput| Finished.");
    }

    public static String getOutputString(OutputLine l) {
        // Convert the HashSet<String> wt_prots to String
        String wt_prots_string = "";
        for (String s : l.wt_prots)
            wt_prots_string = wt_prots_string.concat(s).concat("|");
        // crop last "|"
        if (wt_prots_string.endsWith("|"))
            wt_prots_string = wt_prots_string.substring(0, wt_prots_string.length() - 1);

        // Convert the HashSet<String> sv_prots to String
        String sv_prots_string = "";
        for (String s : l.sv_prots)
            sv_prots_string = sv_prots_string.concat(s).concat("|");
        // crop last "|"
        if (sv_prots_string.endsWith("|"))
            sv_prots_string = sv_prots_string.substring(0, sv_prots_string.length() - 1);

        // Convert the HashSet<String> wt_prots to String
        String wt_string = Integer.toString(l.sv.start + 1);
        for (Intron intr : l.wt)
            wt_string += ":" + Integer.toString(intr.start) + "|" + Integer.toString(intr.end + 1);
        wt_string += ":" + Integer.toString(l.sv.end);

        return l.gene_id + "\t" + l.gene_symbol + "\t" + l.chromosome + "\t" + l.strand + "\t" + l.nprots + "\t" + l.ntrans
                + "\t" + (l.sv.start + 1) + ":" + l.sv.end + "\t" + wt_string + "\t" + wt_prots_string + "\t" + sv_prots_string
                + "\t" + l.min_skipped_exon + "\t" + l.max_skipped_exon + "\t" + (l.min_skipped_base + 1) + "\t" + (l.max_skipped_base + 1);
    }

    public static void checkIfOutputLineIsEmpty(OutputLine l) {
        String errorMessage = "";
        if (l.gene_id.isEmpty())
            errorMessage = errorMessage.concat("gene_id ");
        if (l.gene_symbol.isEmpty())
            errorMessage = errorMessage.concat("gene_symbol ");
        if (l.chromosome.isEmpty())
            errorMessage = errorMessage.concat("chromosome ");
        if (l.strand.isEmpty())
            errorMessage = errorMessage.concat("strand ");
        if (l.nprots == 0)
            errorMessage = errorMessage.concat("nprots ");
        if (l.ntrans == 0)
            errorMessage = errorMessage.concat("ntrans ");
        if (l.sv.start < 1 || l.sv.end < 1)
            errorMessage = errorMessage.concat("sv (start or end) ");
        if (l.wt.isEmpty())
            errorMessage = errorMessage.concat("wt (introns) ");
        if (l.sv_prots.isEmpty())
            errorMessage = errorMessage.concat("sv_prots (protein_ids) ");
        if (l.wt_prots.isEmpty())
            errorMessage = errorMessage.concat("wt_prots (protein_ids) ");
        if (l.min_skipped_exon < 1 || l.min_skipped_exon > 10000000)
            errorMessage = errorMessage.concat("min_skipped_exon ");
        if (l.max_skipped_exon < 1 || l.max_skipped_exon > 10000000)
            errorMessage = errorMessage.concat("min_skipped_exon ");
        if (l.min_skipped_base < 1 || l.min_skipped_base > 10000000)
            errorMessage = errorMessage.concat("min_skipped_base ");
        if (l.max_skipped_base < 1 || l.max_skipped_base > 10000000)
            errorMessage = errorMessage.concat("max_skipped_base ");
        if (errorMessage.length() > 0) {
            System.out.println("|checkIfOutputLineIsEmpty| Empty value in: \'" + errorMessage + "\'. Line: " + getOutputString(l));
        }
    }
}
