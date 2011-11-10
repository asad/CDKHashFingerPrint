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
package helper;

import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Generates fragments from a given molecule.
 * <p/>
 * Fragments are created by splitting the molecule at non-terminal and
 * non-ring bonds.
 * 
 * Updated to recursively process the fragment list.
 *
 * @author Rajarshi Guha <rguha@indiana.edu>
 * @author updated by Syed Asad Rahman<asad@ebi.ac.uk>
 */
public class Fragmenter2 {

    private static int nround = 0;

    /**
     * Split a molecule into fragments.
     * <p/>
     * The method considers bonds as splittable if they are
     * <ul>
     * <li>non-terminal
     * <li>Not in a ring
     * </ul>
     * In addition, the fragments returned will include rings if present.
     * <p/>
     * Make sure to remove hydrogens before calling this method!
     *
     * @param atomContainer The molecule to split
     * @param verbose 
     * @return a list of fragments
     * @throws org.openscience.cdk.exception.CDKException
     *
     */
    public List<IAtomContainer> generateFragments(IAtomContainer atomContainer, boolean verbose)
            throws CDKException {

        nround++;

        ArrayList<IAtomContainer> fragments = new ArrayList<IAtomContainer>();

        if (atomContainer.getBondCount() < 3) {
            return fragments;
        }

        List<IBond> splitableBonds = getSplitableBonds(atomContainer);
        if (splitableBonds.isEmpty()) {
            return fragments;
        }

        if (verbose) {
            System.out.println("Found " + splitableBonds.size() + " bonds to split on");
        }

        for (IBond bond : splitableBonds) {
            List<IAtomContainer> parts = splitMolecule(atomContainer, bond);

            if (fragments.isEmpty()) {
                fragments.addAll(parts);
                continue;
            }

            // make sure we don't add the same fragment twice
            for (IAtomContainer partContainer : parts) {
                if (!fragmentExists(partContainer, fragments)) {
                    fragments.add(partContainer);
                }
            }
        }

        // try and partition the fragments
        List<IAtomContainer> tmp = new ArrayList<IAtomContainer>(fragments);
        for (IAtomContainer fragment : fragments) {
            if (fragment.getBondCount() < 3) {
                continue;
            }
            if (getSplitableBonds(fragment).isEmpty()) {
                continue;
            }

            List<IAtomContainer> frags = generateFragments(fragment, false);
            if (frags.isEmpty()) {
                continue;
            }

            for (IAtomContainer frag : frags) {
                if (frag.getBondCount() < 2) {
                    continue;
                }
                if (!fragmentExists(frag, tmp)) {
                    tmp.add(frag);
                }
            }
        }
        fragments = new ArrayList<IAtomContainer>(tmp);
        return fragments;
    }

    private boolean fragmentExists(IAtomContainer atomContainer, List<IAtomContainer> fragments) {
        boolean present = false;
        for (IAtomContainer f : fragments) {
            if (identicalAtoms(f, atomContainer)) {
                present = true;
                break;
            }
        }
        return present;
    }

    private List<IBond> getSplitableBonds(IAtomContainer atomContainer) throws CDKException {
        // do ring detection
        AllRingsFinder allRingsFinder = new AllRingsFinder();
        IRingSet allRings;
        allRings = allRingsFinder.findAllRings(atomContainer);

        // find the splitable bonds
        ArrayList<IBond> splitableBonds = new ArrayList<IBond>();
        for (IBond bond : atomContainer.bonds()) {

            boolean isInRing = false;
            boolean isTerminal = false;

            // lets see if it's in a ring
            IRingSet rings = allRings.getRings(bond);
            if (rings.getAtomContainerCount() > 0) {
                isInRing = true;
            }

            // lets see if it is a terminal bond
            for (IAtom atom : bond.atoms()) {
                if (atomContainer.getConnectedAtomsCount(atom) == 1) {
                    isTerminal = true;
                    break;
                }
            }

            if (!(isInRing || isTerminal)) {
                splitableBonds.add(bond);
            }
        }
        return splitableBonds;
    }

    private boolean identicalAtoms(IAtomContainer molecule1,
            IAtomContainer molecule2) {
        if (molecule1.getBondCount() != molecule2.getBondCount()
                && molecule1.getAtomCount() != molecule2.getAtomCount()) {
            return false;
        }

        int natom = molecule1.getAtomCount();
        int n = 0;
        for (int i = 0; i < molecule1.getAtomCount(); i++) {
            for (int j = 0; j < molecule2.getAtomCount(); j++) {
                if (molecule1.getAtom(i).equals(molecule2.getAtom(j))) {
                    n++;
                    break;
                }
            }
        }
        return n == molecule1.getAtomCount();
    }

    /*
    in this version we just split a molecule by deleting the bond in question
    If this cases one of the parts to be a single atom, or a single bond we
    ignore that part
     */
    private List<IAtomContainer> xsplitMolecule(IAtomContainer atomContainer, IBond bond) {
        List<IAtomContainer> ret = new ArrayList<IAtomContainer>();

        for (IAtom atom : bond.atoms()) {
            List<IBond> part = new ArrayList<IBond>();
            part.add(bond);
            part = traverse(atomContainer, atom, part);
            part.remove(0); // get rid of the splitting bond

            // this was a terminal bond and we gstarted from the terminal atom, ignore
            if (part.isEmpty()) {
                continue;
            }

            // this was a single bond, ignore
            if (part.size() == 1) {
                continue;
            }

            IAtomContainer partContainer = makeAtomContainer(atom, part);
            ret.add(partContainer);
        }
        return ret;
    }

    // the simplest approach would be to delete the specified bond and then
    // extract the two disconnected fragments. This implies that we have to
    // create an IAtomContainer object each time this method is called. To avoid
    // this we traverse the graph rather than do anything destructive
    private List<IAtomContainer> splitMolecule(IAtomContainer atomContainer,
            IBond bond) {
        List<IAtomContainer> ret = new ArrayList<IAtomContainer>();
        for (IAtom atom: bond.atoms()) {
            List<IBond> part = new ArrayList<IBond>();
            part.add(bond);
            part = traverse(atomContainer, atom, part);

            // at this point we have a partion which contains the bond we
            // split. This partition should actually 2 partitions:
            // - one with the splitting bond
            // - one without the splitting bond
            // note that this will lead to repeated fragments when we  do this
            // with adjacent bonds, so when we gather all the fragments we need
            // to check for repeats
            IAtomContainer partContainer;
            partContainer = makeAtomContainer(atom, part);

            // by checking for more than 2 atoms, we exclude single bond fragments
            // also if a fragment has the same number of atoms as the parent molecule,
            // it is the parent molecule, so we exclude it.
            if (partContainer.getAtomCount() > 2 && partContainer.getAtomCount() != atomContainer.getAtomCount()) {
                ret.add(partContainer);
            }

            part.remove(0);
            partContainer = makeAtomContainer(atom, part);
            if (partContainer.getAtomCount() > 2 && partContainer.getAtomCount() != atomContainer.getAtomCount()) {
                ret.add(partContainer);
            }
        }
        return ret;
    }

    private IAtomContainer makeAtomContainer(IAtom atom, List<IBond> parts) {
        IAtomContainer partContainer = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        partContainer.addAtom(atom);
        for (IBond aBond : parts) {
            for (IAtom bondedAtom : aBond.atoms()) {
                if (!partContainer.contains(bondedAtom)) {
                    partContainer.addAtom(bondedAtom);
                }
            }
            partContainer.addBond(aBond);
        }
        return partContainer;
    }

    private List<IBond> traverse(IAtomContainer atomContainer, IAtom atom,
            List<IBond> bondList) {
        List<IBond> connectedBonds = atomContainer.getConnectedBondsList(atom);
        for (IBond aBond : connectedBonds) {
            if (bondList.contains(aBond)) {
                continue;
            }
            bondList.add(aBond);
            IAtom nextAtom = aBond.getConnectedAtom(atom);
            if (atomContainer.getConnectedAtomsCount(nextAtom) == 1) {
                continue;
            }
            traverse(atomContainer, nextAtom, bondList);
        }
        return bondList;
    }

    public IAtomContainer example4() throws InvalidSmilesException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//		return sp.parseSmiles("c1c(CC2CC2)cc(CNCC)cc1");
        return sp.parseSmiles("c1c(CC2CC2)cc(CNCCC)cc1");
//        return sp.parseSmiles("CCC1CC=C1");
//        return sp.parseSmiles("O=[N+]([O-])C=1C=C(C(O)=C(C=1)Cl)[N+](=O)[O-]");
//        return sp.parseSmiles("SPCN");
    }

    public static void main(String[] args) throws CDKException {
        IAtomContainer mol = null;

        Fragmenter2 fr = new Fragmenter2();

        if (args.length == 0) {
            System.out.println("You can specify a SMILES string on the command line");
            mol = fr.example4();
        }
        if (args.length > 1) {
            System.out.println("Must supply a single SMILES string");
            System.exit(-1);
        } else if (args.length == 1) {
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            mol = sp.parseSmiles(args[0]);
        }
       
        CDKHueckelAromaticityDetector.detectAromaticity(mol);

        List<IAtomContainer> l = fr.generateFragments(mol, true);
        System.out.println("Got " + l.size() + " fragments");
        System.out.println("Needed " + nround + " calls to generateFragments()");
    }
}