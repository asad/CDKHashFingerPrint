/* Copyright (C) 2013       Syed Asad Rahman <S9asad@gmail.com>
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
package com.bioinception.chem.fp.fingerprints.feature;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 * Generates a features of atom for a given {@link IAtomContainer}. Features
 * which are considered are mentioned below.
 *
 * The first seven feature we extract are as follows: (i) number of non-hydrogen
 * connections, (ii) number of non-hydrogen bonds, (iii) atomic numbers, (iv)
 * sign of charge, (v) absolute charge, (vi) number of connected hydrogens, and
 * (vii) atomic numbers of neighboring atoms
 *
 * <p>
 * A features is generated for an AtomContainer with this code:
 * <pre>
 *   Map<IAtom, String> atomInvariants = FeatureGenerator.getAtomInvariants(container);
 * </pre></p>
 *
 * <p>
 * The FeatureGenerator assumes that hydrogens are explicitly given!
 * Furthermore, if pseudo atoms or atoms with malformed symbols are present,
 * their atomic number is taken as one more than the last element currently
 * supported in {@link PeriodicTable}.
 *
 * <p>
 *
 *
 * @cdk.keyword fingerprint
 * @cdk.keyword similarity
 * @cdk.module standard
 * @cdk.githash
 */
public class FeatureGenerator {

    public static void main(String[] args) {
//        String s1 = "N[C@@H](CC1=CC=C(O)C=C1)C(O)=O";
//        String s2 = "OC(=O)C(=O)CC1=CC=C(O)C=C1";

        /*
         Ring-1
         */
        String s1 = "C1CCCCC1";
        /*
         Ring-2
         */
        String s2 = "C1CCC2CCCCC2C1";
//        /*
//         Linear
//         */
//        String s2 = "CCCCCC";
//        /*
//         Rings plus
//         */
//        String s2 = "CC1CCCCC1O";

        final SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        try {
            IAtomContainer parseSmilesQ = sp.parseSmiles(s1);
            Map<IAtom, String> atomInvariants = getAtomInvariants(parseSmilesQ);
            System.out.println("features " + atomInvariants.values());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(FeatureGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }

        try {
            IAtomContainer parseSmilesQ = sp.parseSmiles(s2);
            Map<IAtom, String> atomInvariants = getAtomInvariants(parseSmilesQ);
            System.out.println("features " + atomInvariants.values());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(FeatureGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    /**
     *
     * @param unsortedMap
     * @return
     */
    private static Map<IAtom, String> sortByValue(Map<IAtom, String> unsortedMap) {

        List<Map.Entry<IAtom, String>> list;
        list = new LinkedList<>(unsortedMap.entrySet());

        // Sort list with comparator, to compare the Map values
        Collections.sort(list, new Comparator<Map.Entry<IAtom, String>>() {
            NaturalOrderComparator naturalOrderComparator = new NaturalOrderComparator();

            @Override
            public int compare(Map.Entry<IAtom, String> o1,
                    Map.Entry<IAtom, String> o2) {
                return (naturalOrderComparator.compare(o1.getValue(), o2.getValue()));
            }
        });

        // Convert sorted map back to a Map
        Map<IAtom, String> sortedMap = new LinkedHashMap<>();
        for (Map.Entry<IAtom, String> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }
        return sortedMap;
    }

    /**
     *
     * @param ac
     * @return
     */
    public static Map<IAtom, String> getAtomInvariants(IAtomContainer ac) {
        Map<IAtom, String> map = new HashMap<>();
        for (IAtom a : ac.atoms()) {
            String atomFeature = getAtomFeature(ac, a);
            map.put(a, atomFeature);
        }
        /*
         All the feature are concatenated and
         then rank ordered using the natural ordering for strings
         */
        Map<IAtom, String> sortByValue = sortByValue(map);
        return sortByValue;
    }
    /*
     The first seven feature we extract are as follows: 
     (i) number of non-hydrogen connections, 
     (ii) number of non-hydrogen bonds, 
     (iii) atomic numbers, 
     (iv) sign of charge, 
     (v) absolute charge, 
     (vi) number of connected hydrogens, and 
     (vii) atomic numbers of neighboring atoms
     */

    /**
     *
     * @param ac
     * @param a
     * @return
     */
    public static String getAtomFeature(IAtomContainer ac, IAtom a) {
        List<String> feature = new ArrayList<>();
//        //(0) Symbol of the atom
//        String symbol = a.getSymbol();
//        feature.add("[" + symbol + "]");
        
        //(i) number of non-hydrogen connections
        int nonHConnections = getNonHydrogenConnections(ac.getConnectedAtomsList(a));
        feature.add(nonHConnections + "");
        //System.out.println("feature " + feature);
        //(ii) number of non-hydrogen bonds
        int nonHBond = getNonHydrogenBonds(ac.getConnectedBondsList(a));
        feature.add(nonHBond + "");
        //System.out.println("feature " + feature);
        //(iii) atomic numbers
        int atomicNumber
                = Objects.equals(a.getAtomicNumber(), CDKConstants.UNSET)
                        ? 0 : a.getAtomicNumber();
        int atomNumber = atomicNumber;
        feature.add(atomNumber + "");
        //System.out.println("feature " + feature);
        //(iv) sign of charge
        Double charge
                = Objects.equals(a.getCharge(), CDKConstants.UNSET) ? 0.0 : a.getCharge();
        // Sign of charge
        if (charge < 0) {
            feature.add(1 + "");
        } else {
            feature.add(0 + "");
        }
        //System.out.println("feature " + feature);
        //(v) absolute charge
        int absCharge
                = (int) Math.abs(
                        (Objects.equals(a.getFormalCharge(), CDKConstants.UNSET)
                                ? 0.0 : a.getFormalCharge()));
        feature.add(absCharge + "");
        //System.out.println("feature " + feature);
        //(vi) number of connected hydrogens
        int countHydrogens = getTotalHyderogenCount(ac, a);
        feature.add(countHydrogens + "");
        //(vii) atomic numbers of neighboring atoms
        List<Integer> neighboringAtomicNumbers = getNeighboringAtomicNumbers(ac.getConnectedAtomsList(a));
        for (Integer i : neighboringAtomicNumbers) {
            feature.add(i + "");
        }
        /*
         All the feature are concatenated and 
         then rank ordered using the natural ordering for strings
         */
        Collections.sort(feature, new NaturalOrderComparator());
        StringBuilder sb = new StringBuilder();
        for (String s : feature) {
            sb.append(s);
        }
        return sb.toString();
    }

    /*
     @returns number of connected hydrogens
     */
    private static int getTotalHyderogenCount(IAtomContainer ac, IAtom a) {
        int hydrogenCount = AtomContainerManipulator.countExplicitHydrogens(ac, a);
        if (a.getImplicitHydrogenCount() == null) {
            throw new IllegalArgumentException("an atom had with unknown (null) implicit hydrogens");
        }
        hydrogenCount += a.getImplicitHydrogenCount();
//        System.out.println("hydrogenCount " + hydrogenCount);
        return hydrogenCount;
    }

    /*
     @returns atomic numbers of neighboring atoms
     */
    private static List<Integer> getNeighboringAtomicNumbers(List<IAtom> connectedAtomsList) {
        List<Integer> atomicNumbersList = new ArrayList<>();
        for (IAtom a : connectedAtomsList) {
            if (a.getAtomicNumber() == null) {
                throw new IllegalArgumentException("an atom had with unknown (null) atomic number");
            }
//            System.out.println("getAtomicNumber " + a.getAtomicNumber());
            atomicNumbersList.add(a.getAtomicNumber());
        }
        return atomicNumbersList;
    }

    /*
     @returns number of non-hydrogen connections
     */
    private static int getNonHydrogenConnections(List<IAtom> connectedAtomsList) {
        int i = 0;
        for (IAtom a : connectedAtomsList) {
            if (Elements.HYDROGEN.getSymbol().equals(a.getSymbol())) {
                continue;
            }
            i++;
        }
        return i;
    }
    /*
     @returns number of hydrogen connections
     */

    private static int getNonHydrogenBonds(List<IBond> connectedAtomsList) {
        int i = 0;
        for (IBond b : connectedAtomsList) {
            IAtom a1 = b.getAtom(0);
            IAtom a2 = b.getAtom(1);

            if (Elements.HYDROGEN.getSymbol().equals(a1.getSymbol())
                    || Elements.HYDROGEN.getSymbol().equals(a2.getSymbol())) {
                continue;
            }
            i += b.getOrder().numeric();
        }
        return i;
    }
}
