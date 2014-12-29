/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2011       Syed Asad Rahman <asad@ebi.ac.uk>
 *           
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package benchmark;

import benchmark.helper.Base;
import benchmark.helper.Data;
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
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.algorithm.vflib.substructure.VF2;

/**
 * Test new FP java -jar dist/CDKHashedFingerprint.jar test/data/mol new 2 1000
 * Test CDK default FP java -jar dist/CDKHashedFingerprint.jar test/data/mol cdk
 * 2 1000 Test new FP with ring matcher java -jar dist/CDKHashedFingerprint.jar
 * test/data/mol new 1 1000
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
    private static fingerprints.interfaces.IFingerprinter fingerprint1 = new fingerprints.hashed.HashedFingerprinter(1024);
    private static fingerprints.interfaces.IFingerprinter fingerprint2 = new fingerprints.HashedBloomFingerprinter(1024);

    /**
     * @param args the command line arguments
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws IOException
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException, IOException {
        String s = System.getProperty("user.home") + File.separator + args[0];
        File directory = new File(s);
        int expectedDataSize = 100;
        if (args.length >= 4) {
            expectedDataSize = Integer.valueOf(args[3]);
        }
        System.out.print("\n***************************************\n");

        Map<String, IAtomContainer> molecules
                = readMDLMolecules(directory, expectedDataSize);
        System.out.println("\rTotal number of mols read: " + molecules.size());
        int interval = (int) (0.10 * molecules.size());
        System.out.println("Intervals between data points: " + interval);
        System.out.print("\n***************************************\n");
        System.out.print("\n------------------------------------------------------------------------------\n");
        System.out.print("CASES:" + "\t\t");
        System.out.print("TP:" + "\t");
        System.out.print("FP:" + "\t");
        System.out.print("TN:" + "\t");
        System.out.print("FN:" + "\t");
        System.out.print("ACCURACY:" + "\t");
        /*
         * TRUE POSITIVE RATE
         */
        System.out.print("TPR:" + "\t");
        /*
         * FALSE POSITIVE RATE
         */
        System.out.print("FPR:" + "\t");
        System.out.println("Time (mins): ");
        System.out.print("------------------------------------------------------------------------------\n");

        if (args.length > 2 && args[2].equals("1")) {
            fingerprint1.setRespectRingMatches(true);
            fingerprint2.setRespectRingMatches(true);
        }

        for (int k = 0; k < molecules.size(); k += interval) {
            int counter = 1;
            for (String inchiKey : molecules.keySet()) {
                IAtomContainer ac = molecules.get(inchiKey);
                if (dataMap.containsKey(inchiKey)) {
                    continue;
                }
                try {
                    BitSet hashedFingerPrint = null;
                    if (args.length > 1 && args[1].equals("cdk")) {
                        hashedFingerPrint = getCDKFingerprint(ac).asBitSet();
                        dataMap.put(inchiKey, new Data(hashedFingerPrint, ac));
                    } else if (args.length > 1 && args[1].equals("hash")) {
                        hashedFingerPrint = getHashedFingerprint(ac).asBitSet();
                        dataMap.put(inchiKey, new Data(hashedFingerPrint, ac));
                    } else if (args.length > 1 && args[1].equals("hashbloom")) {
                        hashedFingerPrint = getHashedBloomFingerprint(ac).asBitSet();
                        dataMap.put(inchiKey, new Data(hashedFingerPrint, ac));
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("error in generating fp: " + ac.getID());
                }

                if (counter++ == interval) {
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
                    sub = new VF2(fragment.getAtomContainer(), original.getAtomContainer(), true, true, false);
                    boolean TrueMatch = sub.isSubgraph();

                    if (FPMatch && TrueMatch) {
                        TP++;
                    } else if (FPMatch && !TrueMatch) {
                        FP++;
                    } else if (!FPMatch && TrueMatch) {
                        FN++;
                    } else if (!FPMatch && !TrueMatch) {
                        TN++;
                    }
                    HITS++;
                }
            }

            System.out.print(dataMap.size() + "*" + dataMap.size() + "\t\t");
            System.out.print(TP + "\t");
            System.out.print(FP + "\t");
            System.out.print(TN + "\t");
            System.out.print(FN + "\t");
            System.out.print(getAccuracy() + "\t\t");
            System.out.print(getTPR() + "\t");
            System.out.print(getFPR() + "\t");
            System.out.println(getElapsedTime(startTime));
        }

    }

    private static BigDecimal getFPR() {
        return FP == 0 ? new BigDecimal(0.000) : new BigDecimal(FP).divide(new BigDecimal(FP + TN), 3, BigDecimal.ROUND_HALF_UP);
    }

    private static BigDecimal getTPR() {
        return TP == 0 ? new BigDecimal(0.000) : new BigDecimal(TP).divide(new BigDecimal(TP + FN), 3, BigDecimal.ROUND_HALF_UP);
    }

    private static BigDecimal getAccuracy() {
        return new BigDecimal(TP + TN).divide(new BigDecimal(HITS), 3, BigDecimal.ROUND_HALF_UP);
    }

    private static String getElapsedTime(long startTime) {
        DecimalFormat df = new DecimalFormat("#.##");
        return df.format((System.currentTimeMillis() - startTime) / (1000 * 60.0)); //*60 for hrs
    }

    private static IBitFingerprint getCDKFingerprint(IAtomContainer ac) throws CDKException {
        return cdkFingerprint.getBitFingerprint(ac);
    }

    private static IBitFingerprint getHashedFingerprint(IAtomContainer ac) throws CDKException {
        return fingerprint1.getBitFingerprint(ac);
    }

    private static IBitFingerprint getHashedBloomFingerprint(IAtomContainer ac) throws CDKException {
        return fingerprint2.getBitFingerprint(ac);
    }
}
