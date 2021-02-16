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
import com.bioinception.chem.fp.fingerprints.bi.ScaffoldHashedFingerprinter;
import com.bioinception.chem.fp.fingerprints.cdk.Fingerprinter;
import com.google.common.collect.FluentIterable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.VentoFoggia;

/**
 * Test new FP java -jar dist/CDKHashedFingerprint.jar test/data/mol scaffold
 * 1000
 *
 * Test CDK default FP java -jar dist/CDKHashedFingerprint.jar test/data/mol cdk
 * 1000
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class BenchmarkFingerprint extends Base {

    private static final Map<String, List<Data>> dataMap = new HashMap<String, List<Data>>();
    private static final org.openscience.cdk.fingerprint.IFingerprinter cdkFingerprint = new Fingerprinter(1024);
    private static final org.openscience.cdk.fingerprint.IFingerprinter scaffoldFingerprint = new ScaffoldHashedFingerprinter(1024);

    /**
     * @param args the command line arguments
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws IOException
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException, IOException {
        int expectedDataSize = 100;
        File directory = new File("/Users/asad/github/rhea/mol");
        System.out.println("mol file dir path: " + directory.getAbsolutePath());

        System.out.print("\n***************************************\n");

        Map<String, IAtomContainer> molecules
                = readMDLMolecules(directory, expectedDataSize);
        System.out.println("\rTotal number of mols read: " + molecules.size());
        int interval = (int) (0.10 * molecules.size());
        System.out.println("Intervals between data points: " + interval);
        System.out.print("\n***************************************\n");

        for (int k = 0; k < molecules.size();) {
            int counter = 1;
            k += interval;
            System.out.println("K: " + k + ", counter: " + counter + ", interval: " + interval);

            for (String inchiKey : molecules.keySet()) {
                IAtomContainer ac = molecules.get(inchiKey);
//                System.out.println("counter: " + counter + ", inchiKey " + inchiKey);
                if (!dataMap.containsKey(inchiKey)) {
                    dataMap.put(inchiKey, new ArrayList<>());
                }
                try {
                    BitSet hashedFingerPrint = null;
                    hashedFingerPrint = getCDKFingerprint(ac).asBitSet();
                    Data data = new Data(hashedFingerPrint, ac);
                    data.setFPName("CDK");
                    dataMap.get(inchiKey).add(data);
                    hashedFingerPrint = getScaffoldFingerprint(ac).asBitSet();
                    data = new Data(hashedFingerPrint, ac);
                    data.setFPName("SCAFFOLD");
                    dataMap.get(inchiKey).add(data);
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("error in generating fp: " + ac.getID());
                }
                counter++;
                if (k < counter) {
                    break;
                }
            }

            System.out.println("Data " + dataMap.size());

            long startTime = System.currentTimeMillis();

            Set<String> cdk_scaffold_fp = new TreeSet<>();
            Set<String> cdk_scaffold_graph = new TreeSet<>();
            Set<String> scaffold_cdk_fp = new TreeSet<>();
            Set<String> scaffold_cdk_graph = new TreeSet<>();
            /*
             * Matcher
             */
            BondMatcher bondmatcher = BondMatcher.forOrder();
            AtomMatcher atommatcher = AtomMatcher.forElement();
            for (String key1 : dataMap.keySet()) {
                List<Data> original = dataMap.get(key1);
                for (String key2 : dataMap.keySet()) {
                    if (key1 == null ? key2 == null : key1.equals(key2)) {
                        continue;
                    }
                    List<Data> fragment = dataMap.get(key2);

//                    System.out.println("key1 " + original.get(0).getAtomContainer().getID()
//                            + ", key2:" + fragment.get(1).getAtomContainer().getID());
                    if (original.get(0).getFPName().equals("CDK")
                            && fragment.get(1).getFPName().equals("SCAFFOLD")) {
                        boolean graph_match = FluentIterable.from(
                                VentoFoggia.findSubstructure(original.get(0).getAtomContainer(),
                                        atommatcher, bondmatcher)
                                        .matchAll(fragment.get(1).getAtomContainer())).size() > 0;
                        boolean fpMatch = FingerprinterTool.isSubset(original.get(0).getFingerprint(),
                                fragment.get(1).getFingerprint());

//                        System.out.println(
//                                "key1 " + key1 + ", key2:" + key2
//                                + ", graph_match " + graph_match
//                                + ", fpMatch " + fpMatch);
                        if (graph_match) {
                            cdk_scaffold_graph.add(original.get(0).getAtomContainer().getID() + "_"
                                    + fragment.get(1).getAtomContainer().getID());
                        }

                        if (fpMatch) {
                            cdk_scaffold_fp.add(original.get(0).getAtomContainer().getID() + "_"
                                    + fragment.get(1).getAtomContainer().getID());
                        }
                    } else if (original.get(1).getFPName().equals("SCAFFOLD")
                            && fragment.get(0).getFPName().equals("CDK")) {
                        boolean graph_match = FluentIterable.from(
                                VentoFoggia.findSubstructure(original.get(0).getAtomContainer())
                                        .matchAll(fragment.get(1).getAtomContainer())).size() > 0;
                        boolean fpMatch = FingerprinterTool.isSubset(original.get(0).getFingerprint(),
                                fragment.get(1).getFingerprint());

//                        System.out.println("key1 " + key1 + ", key2:" + key2
//                                + ", graph_match " + graph_match
//                                + ", fpMatch " + fpMatch);
                        if (graph_match) {
                            scaffold_cdk_graph.add(original.get(0).getAtomContainer().getID() + "_"
                                    + fragment.get(1).getAtomContainer().getID());
                        }

                        if (fpMatch) {
                            scaffold_cdk_fp.add(original.get(0).getAtomContainer().getID() + "_"
                                    + fragment.get(1).getAtomContainer().getID());
                        }
                    }

                }
            }

            Set<String> intersection_cdk_graph = new HashSet<>(cdk_scaffold_fp);
            intersection_cdk_graph.retainAll(cdk_scaffold_graph);
            System.out.println("A: CDK and Graph " + intersection_cdk_graph.size());

            Set<String> intersection_scaffold_graph = new HashSet<>(scaffold_cdk_fp);
            intersection_scaffold_graph.retainAll(scaffold_cdk_graph);
            System.out.println("B: Scaffold and Graph " + intersection_scaffold_graph.size());

            Set<String> intersection = new HashSet<>(intersection_cdk_graph);
            intersection.retainAll(intersection_scaffold_graph);
            System.out.println("C: Scaffold & CDK and Graph " + intersection.size());

            System.out.print(dataMap.size() + "*" + dataMap.size() + "\t\t");
            System.out.println(getElapsedTime(startTime));

            /*
             * break the loop
             */
            if (expectedDataSize == dataMap.size()) {
                break;
            }

        }

    }

    private static String getElapsedTime(long startTime) {
        DecimalFormat df = new DecimalFormat("#.##");
        return df.format((System.currentTimeMillis() - startTime) / (1000 * 60.0)); //*60 for hrs
    }

    private static IBitFingerprint getCDKFingerprint(IAtomContainer ac) throws CDKException {
        return cdkFingerprint.getBitFingerprint(ac);
    }

    private static IBitFingerprint getScaffoldFingerprint(IAtomContainer ac) throws CDKException {
        return scaffoldFingerprint.getBitFingerprint(ac);
    }
}
