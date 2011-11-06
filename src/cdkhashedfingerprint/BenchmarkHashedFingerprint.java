/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cdkhashedfingerprint;

import helper.Data;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.algorithm.vflib.substructure.VF2;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class BenchmarkHashedFingerprint extends Base {

    private static long TP;
    private static long FP;
    private static long FN;
    private static long TN;
    private static long HITS;
    private static Map<String, Data> dataMap = new HashMap<String, Data>();
    private static IFingerprinter cdkFingerprint = new org.openscience.cdk.fingerprint.Fingerprinter(1024);
    private static IFingerprinter fingerprint = new fingerprints.Fingerprinter(1024);

    /**
     * @param args the command line arguments
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws IOException  
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException, IOException {
        String s = System.getProperty("user.home") + File.separator + args[0];
        File directory = new File(s);
        int expectedDataSize = 2000;
        Map<String, IAtomContainer> molecules =
                readMDLMolecules(directory, expectedDataSize);
        System.out.println("Total number of mols read: " + molecules.size());
        int interval = (int) (0.10 * molecules.size());
        System.out.println("interval " + interval);
        System.out.print("\n***************************************\n");
        System.out.print("#CASES:" + "\t");
        System.out.print("TP:" + "\t");
        System.out.print("FP:" + "\t");
        System.out.print("TN:" + "\t");
        System.out.print("FN:" + "\t");
        System.out.print("ACCURACY:" + "\t");
        /*TRUE POSITIVE RATE*/
        System.out.print("TPR:" + "\t");
        /*FALSE POSITIVE RATE*/
        System.out.print("FPR:" + "\t");
        System.out.println("Time (mins): ");

        for (int k = interval; k < molecules.size(); k += interval) {
            for (String inchiKey : molecules.keySet()) {
                IAtomContainer ac = molecules.get(inchiKey);
                if (dataMap.containsKey(inchiKey)) {
                    continue;
                }
                try {
                    BitSet hashedFingerPrint;
                    if (args.length > 1 && args[1].equals("cdk")) {
                        hashedFingerPrint = getCDKFingerprint(ac);
                        dataMap.put(inchiKey, new Data(hashedFingerPrint, ac));
                    } else if (args.length > 1 && args[1].equals("new")) {
                        hashedFingerPrint = getNewFingerprint(ac);
                        dataMap.put(inchiKey, new Data(hashedFingerPrint, ac));
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("error in generating fp: " + ac.getID());
                }

                if (k == dataMap.size()) {
                    break;
                }
            }

            TP = 0;
            FP = 0;
            FN = 0;
            TN = 0;
            HITS = 0;

            long startTime = System.currentTimeMillis();
            for (Data fragment : dataMap.values()) {
                for (Data original : dataMap.values()) {
                    boolean FPMatch = FingerprinterTool.isSubset(
                            original.getFingerprint(),
                            fragment.getFingerprint());
                    VF2 sub;
                    if (args[2].equals("1")) {
                        sub = new VF2(true, true);
                    } else {
                        sub = new VF2(true, false);
                    }
                    sub.set(fragment.getAtomContainer(), original.getAtomContainer());
                    boolean TrueMatch = sub.isSubgraph();

                    if (FPMatch && TrueMatch) {
                        TP++;
                    }
                    if (FPMatch && !TrueMatch) {
                        FP++;
                    }
                    if (!FPMatch && TrueMatch) {
                        FN++;
                    }
                    if (!FPMatch && !TrueMatch) {
                        TN++;
                    }
                    HITS++;
                }
            }

            System.out.print(dataMap.size() + "*" + dataMap.size() + "\t");
            System.out.print(TP + "\t");
            System.out.print(FP + "\t");
            System.out.print(TN + "\t");
            System.out.print(FN + "\t");
            System.out.print(getAccuracy() + "\t");
            System.out.print(getTPR() + "\t");
            System.out.print(getFPR() + "\t");
            System.out.println(getElapsedTime(startTime));
        }
        System.out.print("\n***************************************\n");
    }

    private static BigDecimal getFPR() {
        return new BigDecimal(FP).divide(new BigDecimal(FP + TN), 3, BigDecimal.ROUND_HALF_UP);
    }

    private static BigDecimal getTPR() {
        return new BigDecimal(TP).divide(new BigDecimal(TP + FN), 3, BigDecimal.ROUND_HALF_UP);
    }

    private static BigDecimal getAccuracy() {
        return new BigDecimal(TP + TN).divide(new BigDecimal(HITS), 3, BigDecimal.ROUND_HALF_UP);
    }

    private static String getElapsedTime(long startTime) {
        DecimalFormat df = new DecimalFormat("#.##");
        return df.format((System.currentTimeMillis() - startTime) / (1000 * 60.0)); //*60 for hrs
    }

    private static BitSet getCDKFingerprint(IAtomContainer ac) throws CDKException {
        return cdkFingerprint.getFingerprint(ac);
    }

    private static BitSet getNewFingerprint(IAtomContainer ac) throws CDKException {
        return fingerprint.getFingerprint(ac);
    }
}
