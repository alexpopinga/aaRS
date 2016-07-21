package aaRS;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class POAtoFASTA {

    public static void main(String[] args) throws IOException {

        FileReader fileReader = new FileReader(args[0]);
        BufferedReader reader = new BufferedReader(fileReader);

        String line = reader.readLine();
        List<String> seqNames = new ArrayList<String>();

        while (! line.equals("<PRO>")) {
            line = reader.readLine();
        }

        line = reader.readLine();

        while (! line.equals("</PRO>")) {
            String[] s = line.split(" ");
            seqNames.add(s[0]);
            line = reader.readLine();
        }

        StringBuilder[] seqs = new StringBuilder[seqNames.size()];

        for (int i=0; i<seqs.length; i++) {
            seqs[i] = new StringBuilder();
        }

        int alignmentLength = 0;

        while (! line.equals("<NODE>")) {
            line = reader.readLine();
        }

        line = reader.readLine();

        while (! line.equals("</NODE>")) {
            String[] s = line.split(" ");
            int n = s.length / 2;

            for (int i=0; i<n; i++) {
                int seq = Integer.parseInt(s[i*2+1]);
                String aaStr = s[(i+1)*2];
                String aminoAcid = aaStr.split("\\.")[2];

                seqs[seq].append(aminoAcid);
            }

            alignmentLength += 1;

            for (int i=0; i<seqs.length; i++) {
                if (seqs[i].length() < alignmentLength) {
                    seqs[i].append("-");
                }
            }

            line = reader.readLine();
        }

        reader.close();
        PrintWriter writer = new PrintWriter("FASTA.txt", "UTF-8");

        for (int i=0; i<seqs.length; i++) {
            System.out.println(">" + seqNames.get(i));
            System.out.println(seqs[i].toString());

            writer.println(">" + seqNames.get(i));
            writer.println(seqs[i].toString());
        }
        writer.flush();
        writer.close();
    }
}
