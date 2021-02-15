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
package com.bioinception.chem.fp.fingerprints.helper;

import java.io.Serializable;
import java.util.BitSet;
import java.util.Collection;
import org.apache.commons.lang3.builder.HashCodeBuilder;

/**
 * A simple Bloom Filter (see http://en.wikipedia.org/wiki/Bloom_filter)
 *
 * Inspired by the SimpleBloomFilter-class written by Ian Clarke. This
 * implementation provides a more evenly distributed Hash-function by using a
 * proper digest instead of the Java RNG. Many of the changes were proposed in
 * comments in his blog:
 * http://blog.locut.us/2008/01/12/a-decent-stand-alone-java-bloom-filter-implementation/
 *
 * Copyright (c) Ian Clarke Copyright (c) Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @param <T> Bloom filter data/object type.
 */
public class BloomFilter<T> extends RandomNumber implements Serializable {

    private static final long serialVersionUID = 0x2c67671896f80faL;
    protected final int k;
    protected final BitSet bitSet;
    private final int bitSetSize, expectedPatterns;

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final BloomFilter<T> other = (BloomFilter<T>) obj;
        if (this.k != other.k) {
            return false;
        }
        if (this.bitSet != other.bitSet && (this.bitSet == null || !this.bitSet.equals(other.bitSet))) {
            return false;
        }
        if (this.bitSetSize != other.bitSetSize) {
            return false;
        }
        return this.expectedPatterns == other.expectedPatterns;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 17 * hash + this.k;
        hash = 17 * hash + (this.bitSet != null
                ? (new HashCodeBuilder(17, 37).append(this.bitSet).toHashCode()) : 0);
        hash = 17 * hash + this.bitSetSize;
        hash = 17 * hash + this.expectedPatterns;
        return hash;
    }

    /**
     * Calculates a hash code for this class.
     *
     * @return hash code representing the contents of an instance of this class.
     */
    /**
     * Simplified constructor that sets the number of expected elements equal to
     * the number of bits.
     *
     * @param bitArraySize The number of bits in the bit array (often called 'm'
     * in the context of bloom filters).
     */
    public BloomFilter(int bitArraySize) {
        this(bitArraySize, bitArraySize);
    }

    /**
     * You must specify the number of bits in the Bloom Filter, and also you
     * should specify the number of items you expect to add.
     *
     * The latter is used to choose some optimal internal values to minimize the
     * false-positive rate (which can be estimated with
     * expectedFalsePositiveRate()).
     *
     * @param bisetSize The number of bits in the bit array (often called 'm' in
     * the context of bloom filters).
     * @param expectedPatterns The typical number of items you expect to be
     * added to the SimpleBloomFilter (often called 'n').
     */
    public BloomFilter(int bisetSize, int expectedPatterns) {
        this.bitSetSize = bisetSize;
        this.expectedPatterns = expectedPatterns;
        // k = ceil(-log_2(false prob.))
        double falsePositiveProbability = (bisetSize / expectedPatterns);
        this.k = (int) Math.ceil((falsePositiveProbability) * Math.log(2.0));
        bitSet = new BitSet(bisetSize);
    }

    /**
     * Calculate the probability of a false positive given the specified number
     * of inserted elements.
     *
     * Calculates the approximate probability of the contains() method returning
     * true for an object that had not previously been inserted into the bloom
     * filter. This is known as the "false positive probability".
     *
     * @return The estimated false positive rate
     */
    public double getExpectedFalsePositiveProbability() {
        // (1 - e^(-k * n / m)) ^ k
        return Math.pow((1 - Math.exp(-k * (double) expectedPatterns / (double) bitSetSize)), k);
    }

    /**
     * Returns the value chosen for K.<br /> <br /> K is the optimal number of
     * hash functions based on the size of the Bloom filter and the expected
     * number of inserted elements.
     *
     * @return optimal k.
     */
    public int getK() {
        return k;
    }

    /**
     * Sets all bits to false in the Bloom filter.
     */
    public void clear() {
        bitSet.clear();
    }


    /*
     * @return This method will always return false
     *
     * @see java.util.Set#add(java.lang.Object)
     */
    public boolean add(T o) {
        int toHashCode = new HashCodeBuilder(17, 37).append(o).toHashCode();
        for (int i = 0; i < k; i++) {
            int position = (int) generateMersenneTwisterRandomNumber(bitSetSize, toHashCode);
            setBit(position, true);
        }
        return false;
    }

    /**
     * @param c
     * @return This method will always return false
     */
    public boolean addAll(Collection<? extends T> c) {
        c.forEach(t -> {
            add(t);
        });
        return false;
    }

    /**
     * Returns true if the element could have been inserted into the Bloom
     * filter. Use getFalsePositiveProbability() to calculate the probability of
     * this being correct.
     *
     * @param o object to check.
     * @return False indicates that t was definitely not added to this Bloom
     * Filter, true indicates that it probably was. The probability can be
     * estimated using the getExpectedFalsePositiveProbability() method.
     */
    public boolean contains(Object o) {
        int toHashCode = new HashCodeBuilder(17, 37).append(o).toHashCode();
        for (int x = 0; x < k; x++) {
            int position = (int) generateMersenneTwisterRandomNumber(bitSetSize, toHashCode);
            if (!getBit(position)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns true if all the elements of a Collection could have been inserted
     * into the Bloom filter. Use getFalsePositiveProbability() to calculate the
     * probability of this being correct.
     *
     * @param c elements to check.
     * @return true if all the elements in c could have been inserted into the
     * Bloom filter.
     */
    public boolean containsAll(Collection<? extends T> c) {
        return c.stream().noneMatch(o -> (!contains(o)));
    }

    public BitSet toBitSet() {
        BitSet result = new BitSet(bitSetSize);
        result.or(bitSet);
        return result;
    }

    public int getBitArraySize() {
        return bitSetSize;
    }

    public Object[] toArray() {
        int size = getBitArraySize();
        Object a[] = new Object[size];
        int id = 0;
        for (int i = 0; i < size; i++) {
            a[i] = bitSet.get(i) ? 1 : 0;
        }
        return a;
    }

    /**
     * Read a single bit from the Bloom filter.
     *
     * @param bit the bit to read.
     * @return true if the bit is set, false if it is not.
     */
    public boolean getBit(int bit) {
        return bitSet.get(bit);
    }

    /**
     * Set a single bit in the Bloom filter.
     *
     * @param bit is the bit to set.
     * @param value If true, the bit is set. If false, the bit is cleared.
     */
    public void setBit(int bit, boolean value) {
        bitSet.set(bit, value);
    }

    /**
     *
     * @param other
     * @return
     */
    public int intersect(BloomFilter<T> other) {
        BitSet intersection = (BitSet) this.bitSet.clone();
        intersection.and(other.bitSet);
        return (intersection.cardinality() / bitSetSize);
    }

    public BitSetIterable IterableBitSet() {
        return new BitSetIterable(bitSet);
    }
}
