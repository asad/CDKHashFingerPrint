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
package com.bioinception.chem.fp.benchmark.helper;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sf.jniinchi.INCHI_RET;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Base {

    public Base() {
    }

    /**
     *
     * @param ac
     * @return
     * @throws CDKException
     */
    public static String generateInchiKey(IAtomContainer ac) throws CDKException {
        InChIGenerator gen = InChIGeneratorFactory.getInstance().getInChIGenerator(ac);
        if (gen.getReturnStatus() != INCHI_RET.OKAY) {
            //System.err.println("inchi failed: " + gen.getMessage());
            return null;
        }
        return gen.getInchiKey();
    }

    /**
     *
     * @param dir
     * @param cutoff
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws IOException
     */
    public static Map<String, IAtomContainer> readMDLMolecules(File dir, int cutoff) throws FileNotFoundException, CDKException, IOException {
        System.out.println("\nReading Files: ");
        Map<String, IAtomContainer> inchiMolMap = new HashMap<String, IAtomContainer>();
        if (dir.isDirectory()) {
            File[] listFiles = dir.listFiles();
            for (File fileIndex : listFiles) {
                if (fileIndex.isFile() && fileIndex.getName().contains(".mol")) {
                    MDLV2000Reader reader = new MDLV2000Reader(new FileReader(fileIndex));
                    IAtomContainer ac = (IAtomContainer) reader.read(new AtomContainer());
                    try {
                        initializeMolecule(ac);
                        ac.setID((fileIndex.getName().split(".mol"))[0]);
                        String inchiKey = generateInchiKey(ac);
                        if (inchiKey == null || ac.getAtomCount() < 3) {
                            continue;
                        }
                        inchiMolMap.put(inchiKey, ac);
                    } catch (Exception ex) {
                        Logger.getLogger(Base.class.getName()).log(Level.SEVERE, null, ex);
                    }

                    if ((cutoff - inchiMolMap.size()) == 0) {
                        break;
                    } else {
                        System.out.print("\r# " + (cutoff - inchiMolMap.size()));
                    }
                    reader.close();
                }
            }
        }
        return inchiMolMap;
    }

    /**
     * Prepare the target molecule for analysis.
     * <p/>
     * We perform ring perception and aromaticity detection and set up the
     * appropriate properties. Right now, this function is called each time we
     * need to do a query and this is inefficient.
     *
     * @throws CDKException if there is a problem in ring perception or
     * aromaticity detection, which is usually related to a timeout in the ring
     * finding code.
     */
    private static Integer initializeMolecule(IAtomContainer atomContainer) throws CDKException {
        Integer hashRings = 0;
        Map<String, Integer> valencesTable = new HashMap<String, Integer>();
        valencesTable.put("H", 1);
        valencesTable.put("Li", 1);
        valencesTable.put("Be", 2);
        valencesTable.put("B", 3);
        valencesTable.put("C", 4);
        valencesTable.put("N", 5);
        valencesTable.put("O", 6);
        valencesTable.put("F", 7);
        valencesTable.put("Na", 1);
        valencesTable.put("Mg", 2);
        valencesTable.put("Al", 3);
        valencesTable.put("Si", 4);
        valencesTable.put("P", 5);
        valencesTable.put("S", 6);
        valencesTable.put("Cl", 7);
        valencesTable.put("K", 1);
        valencesTable.put("Ca", 2);
        valencesTable.put("Ga", 3);
        valencesTable.put("Ge", 4);
        valencesTable.put("As", 5);
        valencesTable.put("Se", 6);
        valencesTable.put("Br", 7);
        valencesTable.put("Rb", 1);
        valencesTable.put("Sr", 2);
        valencesTable.put("In", 3);
        valencesTable.put("Sn", 4);
        valencesTable.put("Sb", 5);
        valencesTable.put("Te", 6);
        valencesTable.put("I", 7);
        valencesTable.put("Cs", 1);
        valencesTable.put("Ba", 2);
        valencesTable.put("Tl", 3);
        valencesTable.put("Pb", 4);
        valencesTable.put("Bi", 5);
        valencesTable.put("Po", 6);
        valencesTable.put("At", 7);
        valencesTable.put("Fr", 1);
        valencesTable.put("Ra", 2);
        valencesTable.put("Cu", 2);
        valencesTable.put("Mn", 2);
        valencesTable.put("Co", 2);

        // do all ring perception
        AllRingsFinder ringFinder = null;
        if (ringFinder == null) {
            ringFinder = new AllRingsFinder();
        }

        IRingSet allRings;
        try {
            allRings = ringFinder.findAllRings(atomContainer);
            hashRings = allRings.getAtomContainerCount() == 0 ? 0 : 1;
        } catch (CDKException e) {
            throw new CDKException(e.toString(), e);
        }

        // sets SSSR information
        SSSRFinder finder = new SSSRFinder(atomContainer);
        IRingSet sssr = finder.findEssentialRings();

        for (IAtom atom : atomContainer.atoms()) {

            // add a property to each ring atom that will be an array of
            // Integers, indicating what size ring the given atom belongs to
            // Add SSSR ring counts
            if (allRings.contains(atom)) { // it's in a ring
                atom.setFlag(CDKConstants.ISINRING, true);
                atom.setFlag(CDKConstants.ISALIPHATIC, false);
                // lets find which ring sets it is a part of
                List<Integer> ringsizes = new ArrayList<Integer>();
                IRingSet currentRings = allRings.getRings(atom);
                int min = 0;
                for (int i = 0; i < currentRings.getAtomContainerCount(); i++) {
                    int ringSize = currentRings.getAtomContainer(i).getAtomCount();
                    if (min > ringSize) {
                        min = ringSize;
                    }
                    ringsizes.add(ringSize);
                }
                atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
                atom.setProperty(CDKConstants.SMALLEST_RINGS, sssr.getRings(atom));
            } else {
                atom.setFlag(CDKConstants.ISINRING, false);
                atom.setFlag(CDKConstants.ISALIPHATIC, true);
            }

            // determine how many rings bonds each atom is a part of
            int hCount;
            if (atom.getImplicitHydrogenCount() == CDKConstants.UNSET) {
                hCount = 0;
            } else {
                hCount = atom.getImplicitHydrogenCount();
            }

            List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);
            int total = hCount + connectedAtoms.size();
            for (IAtom connectedAtom : connectedAtoms) {
                if (connectedAtom.getSymbol().equals("H")) {
                    hCount++;
                }
            }
            atom.setProperty(CDKConstants.TOTAL_CONNECTIONS, total);
            atom.setProperty(CDKConstants.TOTAL_H_COUNT, hCount);

            if (valencesTable.get(atom.getSymbol()) != null) {
                int formalCharge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
                atom.setValency(valencesTable.get(atom.getSymbol()) - formalCharge);
            }
        }

        for (IBond bond : atomContainer.bonds()) {
            if (allRings.getRings(bond).getAtomContainerCount() > 0) {
                bond.setFlag(CDKConstants.ISINRING, true);
                bond.setFlag(CDKConstants.ISALIPHATIC, false);
            }
        }

        for (IAtom atom : atomContainer.atoms()) {
            List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);

            int counter = 0;
            IAtom any;
            for (IAtom connectedAtom : connectedAtoms) {
                any = connectedAtom;
                if (any.getFlag(CDKConstants.ISINRING)) {
                    counter++;
                }
            }
            atom.setProperty(CDKConstants.RING_CONNECTIONS, counter);
        }

        // check for atomaticity
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
            Aromaticity.cdkLegacy().apply(atomContainer);
        } catch (CDKException e) {
            throw new CDKException(e.toString(), e);
        }
        return hashRings;
    }
}
