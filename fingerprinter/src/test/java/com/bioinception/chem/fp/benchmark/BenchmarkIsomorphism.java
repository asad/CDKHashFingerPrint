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
package com.bioinception.chem.fp.benchmark;

import com.bioinception.chem.fp.benchmark.helper.Base;
import static com.bioinception.chem.fp.benchmark.helper.Base.readMDLMolecules;
import com.bioinception.chem.fp.benchmark.helper.Data;
import com.google.common.collect.FluentIterable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.VentoFoggia;

/**
 * Test new FP java -jar dist/CDKHashedFingerprint.jar test/data/mol new 2 1000
 * Test CDK default FP java -jar dist/CDKHashedFingerprint.jar test/data/mol cdk
 * 2 1000 Test new FP with ring matcher java -jar dist/CDKHashedFingerprint.jar
 * test/data/mol new 1 1000
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class BenchmarkIsomorphism extends Base {

    private static long TP;
    private static long FP;
    private static long FN;
    private static long TN;
    private static long HITS;
    private static final Map<String, Data> dataMap = new HashMap<String, Data>();

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
        if (args.length >= 2) {
            expectedDataSize = Integer.valueOf(args[1]);
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

        for (int k = 0; k < molecules.size(); k += interval) {
            int counter = 1;
            for (String inchiKey : molecules.keySet()) {
                IAtomContainer ac = molecules.get(inchiKey);
                if (dataMap.containsKey(inchiKey)) {
                    continue;
                }
                try {
                    dataMap.put(inchiKey, new Data(ac));
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("error in generating graph: " + ac.getID());
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
            BondMatcher bondmatcher = BondMatcher.forOrder();
            AtomMatcher atommatcher = AtomMatcher.forAny();

            long startTime = System.currentTimeMillis();
            dataMap.values().forEach(fragment -> {
                dataMap.values().stream().map(original -> {
                    int count_bond_match = FluentIterable.from(
                            VentoFoggia.findIdentical(original.getAtomContainer(),
                                    atommatcher, bondmatcher)
                                    .matchAll(fragment.getAtomContainer())).size();
                    int count_match = FluentIterable.from(
                            VentoFoggia.findIdentical(original.getAtomContainer())
                                    .matchAll(fragment.getAtomContainer())).size();
                    boolean trueMatch = count_bond_match > 0;
                    boolean bondMatch = count_match > 0;
                    if (trueMatch && bondMatch) {
                        TP++;
                    } else if (trueMatch && !bondMatch) {
//                        System.out.println("\nQ " + fragment.getAtomContainer().getID());
//                        System.out.println("T " + original.getAtomContainer().getID());
                        FP++;
                    } else if (!trueMatch && bondMatch) {
                        FN++;
                    } else if (!trueMatch && !bondMatch) {
                        TN++;
                    }
                    return original;
                }).forEachOrdered(_item -> {
                    HITS++;
                });
            });

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
}